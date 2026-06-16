#!/bin/bash --login
# =============================================================================
# Experiment 5 — Ablation Study   (Q4: Does the system design matter?)
# =============================================================================
# Roadmap §3.5.7 (Exp 5). Each design decision is toggled by a compile-time or
# env flag and measured in isolation.
#
# WHY THIS REWRITE RUNS AT HPRC SCALE
# -----------------------------------
# The Phase-2/3 ablation at n=10K returned NULL results for A1 (vEB vs BFS) and
# A2/A3 (steal-largest vs steal-random): at 10K the working set (~56 MB) fits in
# the EPYC 7763's 256 MB L3, so cache layout is irrelevant, and per-task DP work
# (k^3, k~2-3) is below scheduler overhead, so steal policy is irrelevant. Those
# nulls are EXPECTED and interpretable — but a BIBM paper needs the effects to
# actually appear. They only appear once the snarl-tree working set exceeds L3,
# i.e. at HPRC scale (~10M snarls). This script therefore runs every ablation at
# BOTH n=10K (control, shows the null) AND n=10M (treatment, shows the effect).
#
# THE FIVE ABLATIONS
#   A1  vEB layout vs BFS layout           -> L3 miss rate (perf) + wall
#   A2  steal-largest vs steal-random      -> wall on standard graph
#   A3  steal-largest vs steal-random      -> wall on heavy-tailed graph (seed 999)
#   A4  MinimizerIndex vs naive O(|V|^2)   -> anchor-detection wall + recall
#   A5  Hungarian(LAPJV) vs brute-force     -> us/op at k in {2,3,5,7,10,20,50}
#
# Each ablation toggles via the documented flags already in the codebase:
#   -DG2G_DISABLE_VEB=1              (A1 BFS layout)
#   -DG2G_STEAL_RANDOM=1            (A2/A3 random steal)
#   -DG2G_DISABLE_MINIMIZER_INDEX=1 (A4 naive anchor)
#
# OUTPUTS
#   ablation_results.tsv   all of A1-A4 (ablation,variant,n,run,wall,l3_miss,rss)
#   matching_bench.tsv     A5 (k, us_per_op, run)
#   ablation_summary.tsv   mean +/- stderr per (ablation, variant, n)
#
# SUBMIT:  sbatch slurm/scripts/exp5_ablation.sh
# =============================================================================

#SBATCH --job-name=g2g_exp5
#SBATCH --output=slurm/logs/exp5-%j.out
#SBATCH --error=slurm/logs/exp5-%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=512GB
#SBATCH --constraint=cpu_amd_epyc_7763

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${G2G_REPO:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
BUILD_BASE="$REPO/build"
GEN="$REPO/slurm/synth/gen_hprc_scale_gfa.py"
OUTDIR="$REPO/results/exp5"
GRAPHDIR="${G2G_GRAPHDIR:-$REPO/results/exp3/graphs}"
RUNS="${G2G_RUNS:-5}"
NCPU="$(nproc)"

# Control (fits in L3) + treatment (exceeds L3). Drop 10M if RAM-constrained.
N_CONTROL="${G2G_N_CONTROL:-10000}"
N_TREATMENT="${G2G_N_TREATMENT:-10000000}"
# Size at which the MinimizerIndex path is triggered (|VA|*|VB| > 1e6).
N_ANCHOR="${G2G_N_ANCHOR:-200000}"
# Heavy-tailed graph for A3 (more skewed allele distribution).
N_HEAVY="${G2G_N_HEAVY:-10000000}"
HEAVY_ALPHA="${G2G_HEAVY_ALPHA:-1.0}"   # smaller alpha => heavier tail

mkdir -p "$OUTDIR" "$GRAPHDIR" "$REPO/slurm/logs"
source "$SCRIPT_DIR/lib_measure.sh"

if [[ -n "${CONDA_BASE:-}" && -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
fi
conda activate "${G2G_ENV:-phase5}" 2>/dev/null || log "WARN: conda env not activated"

capture_env "$OUTDIR/env.txt"
log "exp5 start  control_n=$N_CONTROL  treatment_n=$N_TREATMENT  perf=$HAVE_PERF"

RESULTS="$OUTDIR/ablation_results.tsv"
BENCH="$OUTDIR/matching_bench.tsv"
tsv_init "$RESULTS" "ablation\tvariant\tn_snarls\trun\twall_sec\tpeak_rss_kb\tl3_miss_rate"
tsv_init "$BENCH"   "k\tus_per_op\trun"

# ── Build a g2g variant with optional extra flags, return its path ───────────
build_variant() {
    local suffix="$1"; local extra="${2:-}"
    local bdir="$BUILD_BASE/$suffix"
    if command -v cmake >/dev/null 2>&1; then
        cmake -B "$bdir" -DCMAKE_BUILD_TYPE=Release \
              -DCMAKE_CXX_COMPILER="$(which g++)" -DCMAKE_C_COMPILER="$(which gcc)" \
              -DCMAKE_CXX_FLAGS="-O3 -march=native ${extra}" "$REPO" >/dev/null 2>&1
        cmake --build "$bdir" -j"$NCPU" --target g2g >/dev/null 2>&1
    fi
    echo "$bdir/g2g"
}

# ── Ensure a graph pair exists ───────────────────────────────────────────────
ensure_graph() {
    local n="$1"; local seed="${2:-42}"; local alpha="${3:-2.0}"; local tag="${4:-}"
    local ga gb
    if [[ -n "$tag" ]]; then ga="$GRAPHDIR/ga_${tag}.gfa"; gb="$GRAPHDIR/gb_${tag}.gfa"
    else ga="$GRAPHDIR/ga_${n}.gfa"; gb="$GRAPHDIR/gb_${n}.gfa"; fi
    if [[ ! -s "$ga" || ! -s "$gb" ]]; then
        log "  generating n=$n seed=$seed alpha=$alpha tag='$tag'..."
        python3 "$GEN" --mode size --n_snarls "$n" --out_dir "$GRAPHDIR" \
            --seed "$seed" --alpha "$alpha"
        if [[ -n "$tag" ]]; then
            mv "$GRAPHDIR/ga_${n}.gfa" "$ga"; mv "$GRAPHDIR/gb_${n}.gfa" "$gb"
            mv "$GRAPHDIR/truth_${n}.tsv" "$GRAPHDIR/truth_${tag}.tsv" 2>/dev/null || true
        fi
    fi
    echo "$ga $gb"
}

# =============================================================================
# A1 — vEB layout vs BFS layout  (control + treatment)
# =============================================================================
log "=== A1: vEB vs BFS layout ==="
G2G_VEB=$(build_variant "veb")
G2G_BFS=$(build_variant "bfs" "-DG2G_DISABLE_VEB=1")
[[ -x "$G2G_VEB" && -x "$G2G_BFS" ]] || { log "FATAL: A1 variants missing"; }

for N in $N_CONTROL $N_TREATMENT; do
    read -r GA GB <<< "$(ensure_graph "$N")"
    for RUN in $(seq 1 "$RUNS"); do
        for V in veb bfs; do
            [[ "$V" == "veb" ]] && BIN="$G2G_VEB" || BIN="$G2G_BFS"
            # Sequential to isolate the layout effect from parallelism.
            m=$(measure_full "$BIN" "$GA" "$GB" --sequential --gaf /dev/null)
            read -r wall peak _ _ _ ec <<< "$m"
            L3=$(l3_miss_rate "$BIN" "$GA" "$GB" --sequential --gaf /dev/null)
            printf 'A1\t%s\t%s\t%s\t%s\t%s\t%s\n' "$V" "$N" "$RUN" "$wall" "$peak" "$L3" >> "$RESULTS"
            log "  A1 $V n=$N run=$RUN wall=${wall}s rss=${peak}kB l3=$L3"
        done
    done
done

# =============================================================================
# A2/A3 — steal-largest vs steal-random  (standard + heavy-tailed)
# =============================================================================
log "=== A2/A3: steal-largest vs steal-random ==="
G2G_SL=$(build_variant "steal_largest")
G2G_SR=$(build_variant "steal_random" "-DG2G_STEAL_RANDOM=1")

# A2: standard graph at control + treatment scale
for N in $N_CONTROL $N_TREATMENT; do
    read -r GA GB <<< "$(ensure_graph "$N")"
    for RUN in $(seq 1 "$RUNS"); do
        for V in steal_largest steal_random; do
            [[ "$V" == "steal_largest" ]] && BIN="$G2G_SL" || BIN="$G2G_SR"
            m=$(measure_full "$BIN" "$GA" "$GB" --threads "$NCPU" --gaf /dev/null)
            read -r wall peak _ _ _ ec <<< "$m"
            printf 'A2\t%s\t%s\t%s\t%s\t%s\tNA\n' "$V" "$N" "$RUN" "$wall" "$peak" >> "$RESULTS"
            log "  A2 $V n=$N run=$RUN wall=${wall}s"
        done
    done
done

# A3: heavy-tailed graph (alpha=1.0) at treatment scale — where steal-largest
# should finally beat steal-random because a few snarls dominate the work.
log "  A3: heavy-tailed graph (alpha=$HEAVY_ALPHA)..."
read -r GHA GHB <<< "$(ensure_graph "$N_HEAVY" 999 "$HEAVY_ALPHA" "heavy_${N_HEAVY}")"
for RUN in $(seq 1 "$RUNS"); do
    for V in steal_largest steal_random; do
        [[ "$V" == "steal_largest" ]] && BIN="$G2G_SL" || BIN="$G2G_SR"
        m=$(measure_full "$BIN" "$GHA" "$GHB" --threads "$NCPU" --gaf /dev/null)
        read -r wall peak _ _ _ ec <<< "$m"
        printf 'A3\t%s\theavy_%s\t%s\t%s\t%s\tNA\n' "$V" "$N_HEAVY" "$RUN" "$wall" "$peak" >> "$RESULTS"
        log "  A3 $V heavy n=$N_HEAVY run=$RUN wall=${wall}s"
    done
done

# =============================================================================
# A4 — MinimizerIndex vs naive anchor detection
# =============================================================================
log "=== A4: MinimizerIndex vs naive ==="
G2G_IDX=$(build_variant "minimizer_idx")
G2G_NAIVE=$(build_variant "naive_anchor" "-DG2G_DISABLE_MINIMIZER_INDEX=1")
read -r GAA GAB <<< "$(ensure_graph "$N_ANCHOR")"
for RUN in $(seq 1 "$RUNS"); do
    for V in minimizer_idx naive_anchor; do
        [[ "$V" == "minimizer_idx" ]] && BIN="$G2G_IDX" || BIN="$G2G_NAIVE"
        m=$(measure_full "$BIN" "$GAA" "$GAB" --sequential --gaf /dev/null)
        read -r wall peak _ _ _ ec <<< "$m"
        printf 'A4\t%s\t%s\t%s\t%s\t%s\tNA\n' "$V" "$N_ANCHOR" "$RUN" "$wall" "$peak" >> "$RESULTS"
        log "  A4 $V n=$N_ANCHOR run=$RUN wall=${wall}s"
    done
done

# =============================================================================
# A5 — Hungarian (LAPJV) vs brute-force matching micro-benchmark
# =============================================================================
log "=== A5: matching micro-benchmark ==="
BENCH_SRC="$REPO/slurm/synth/matching_bench.cpp"
BENCH_BIN="$OUTDIR/matching_bench"

# The bench source is shipped in slurm/synth/. If absent, write it.
if [[ ! -f "$BENCH_SRC" ]]; then
    log "  matching_bench.cpp not found; writing canonical version"
    cat > "$BENCH_SRC" <<'CPPSRC'
/**
 * matching_bench.cpp — micro-benchmark for allele_match dispatch (Exp 5, A5).
 * Usage: matching_bench <k> <n_trials> [brute]
 *   - default: times allele_match() (Hungarian/LAPJV for k<=30, greedy above)
 *   - "brute": forces brute_match() (valid only k<=8) for the crossover figure
 * Prints:  <k>\t<microseconds_per_op>
 */
#include <g2g/matching.hpp>
#include <g2g/sequence.hpp>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace g2g;
using Clock = std::chrono::high_resolution_clock;

static std::string rand_dna(int len, std::mt19937& rng) {
    static const char B[] = "ACGT";
    std::string s(len, 'A');
    std::uniform_int_distribution<int> d(0, 3);
    for (char& c : s) c = B[d(rng)];
    return s;
}

int main(int argc, char** argv) {
    if (argc < 3) { std::cerr << "Usage: matching_bench <k> <n_trials> [brute]\n"; return 1; }
    int k = std::atoi(argv[1]);
    int n_trials = std::atoi(argv[2]);
    bool force_brute = (argc > 3 && std::strcmp(argv[3], "brute") == 0);

    std::mt19937 rng(42);
    MatchParams params;

    std::vector<std::string> sa(k), sb(k);
    for (auto& s : sa) s = rand_dna(50, rng);
    for (auto& s : sb) s = rand_dna(50, rng);
    std::vector<std::string_view> va(sa.begin(), sa.end()), vb(sb.begin(), sb.end());

    auto run_one = [&]() {
        if (force_brute) return brute_match(va, vb, params);
        return allele_match(va, vb, params);
    };

    for (int t = 0; t < 10; ++t) run_one();   // warm up

    auto t0 = Clock::now();
    for (int t = 0; t < n_trials; ++t) run_one();
    auto t1 = Clock::now();

    double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
    std::cout << k << "\t" << (us / n_trials) << "\n";
    return 0;
}
CPPSRC
fi

# Compile the bench against the steal_largest static lib (any variant works).
LIBG2G=""
for cand in "$BUILD_BASE/steal_largest/libg2g_lib.a" "$BUILD_BASE/libg2g_lib.a" \
            "$BUILD_BASE/veb/libg2g_lib.a"; do
    [[ -f "$cand" ]] && { LIBG2G="$cand"; break; }
done
if [[ -n "$LIBG2G" ]] && command -v g++ >/dev/null 2>&1; then
    g++ -std=c++20 -O3 -march=native -I"$REPO/include" \
        "$BENCH_SRC" "$LIBG2G" -pthread -fopenmp -o "$BENCH_BIN" 2>/dev/null \
        && log "  matching_bench compiled" || log "  WARN: bench compile failed"
else
    log "  WARN: no libg2g_lib.a found; skipping A5"
fi

if [[ -x "$BENCH_BIN" ]]; then
    # allele_match (LAPJV/greedy) across the full k range
    for K in 2 3 5 7 10 20 50; do
        for RUN in $(seq 1 "$RUNS"); do
            LINE=$("$BENCH_BIN" "$K" 100000)
            printf '%s\t%s\n' "$LINE" "$RUN" >> "$BENCH"
            log "  A5 lapjv k=$K run=$RUN: ${LINE} us/op"
        done
    done
    # brute-force crossover (only up to k=8 where it is feasible)
    BENCH_BRUTE="$OUTDIR/matching_bench_brute.tsv"
    tsv_init "$BENCH_BRUTE" "k\tus_per_op\trun"
    for K in 2 3 5 7 8; do
        for RUN in $(seq 1 "$RUNS"); do
            # Fewer trials at high k since brute is slow.
            TRIALS=$(( K <= 5 ? 10000 : 500 ))
            LINE=$("$BENCH_BIN" "$K" "$TRIALS" brute)
            printf '%s\t%s\n' "$LINE" "$RUN" >> "$BENCH_BRUTE"
            log "  A5 brute k=$K run=$RUN: ${LINE} us/op"
        done
    done
fi

# =============================================================================
# SUMMARY
# =============================================================================
log "computing ablation summary..."
python3 - "$RESULTS" "$OUTDIR/ablation_summary.tsv" <<'PYEOF'
import sys, csv, math
from collections import defaultdict
infile, outfile = sys.argv[1], sys.argv[2]
rows = defaultdict(lambda: {"wall":[], "l3":[]})
with open(infile) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        key=(r["ablation"], r["variant"], r["n_snarls"])
        try: rows[key]["wall"].append(float(r["wall_sec"]))
        except: pass
        if r.get("l3_miss_rate","NA") not in ("NA",""):
            try: rows[key]["l3"].append(float(r["l3_miss_rate"]))
            except: pass
def mean(v): return sum(v)/len(v) if v else 0.0
def se(v):
    if len(v)<2: return 0.0
    m=mean(v); return math.sqrt(sum((x-m)**2 for x in v)/(len(v)-1)/len(v))
with open(outfile,"w") as f:
    f.write("ablation\tvariant\tn_snarls\tmean_wall_sec\tstderr_wall_sec\tmean_l3_miss\tn_runs\n")
    for key in sorted(rows):
        w=rows[key]["wall"]; l3=rows[key]["l3"]
        l3s = f"{mean(l3):.4f}" if l3 else "NA"
        f.write("\t".join(list(key)+[f"{mean(w):.4f}",f"{se(w):.4f}",l3s,str(len(w))])+"\n")
print(f"[exp5] wrote {outfile}")
PYEOF

log "=== Exp5 DONE ==="
log "  ablations:  $RESULTS"
log "  matching:   $BENCH"
log "  summary:    $OUTDIR/ablation_summary.tsv"
