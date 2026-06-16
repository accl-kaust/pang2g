#!/bin/bash --login
# =============================================================================
# Experiment 4 — Parallel Scaling   (Q2)
# =============================================================================
# Roadmap §3.5.6 (Exp 4). Produces the speedup curve + parallel efficiency that
# every parallel-systems paper is expected to show (Figure 3).
#
# THE HONEST-SCALING PROBLEM THIS SCRIPT IS BUILT TO HANDLE
# --------------------------------------------------------
# The Phase-2 100K run showed 1->2 threads gave only 1.002x speedup. Diagnosis:
# at 100K snarls the parallel DP is NOT the bottleneck — GFA parsing, Minimizer
# index construction, and snarl decomposition (all SEQUENTIAL) dominate. Amdahl's
# law then caps speedup at ~1x because <0.2% of runtime is parallelisable.
#
# This script confronts that directly by measuring TWO speedup curves:
#   (A) END-TO-END speedup  = T1_total / Tp_total   (what a user sees)
#   (B) DP-PHASE   speedup  = T1_dp    / Tp_dp       (the runtime's true scaling)
# The gap between (A) and (B) IS the Amdahl story, and reporting both is what a
# reviewer needs to accept the parallel contribution. We push graph size up to
# where DP dominates (>=1M snarls) so (A) and (B) converge.
#
# Per thread count we record, from --verbose-timing and RunStats:
#   wall_total, phase_io, phase_index, phase_decompose, phase_dp,
#   tasks_executed (MUST be constant across p — scheduler correctness check),
#   peak_rss, cpu_pct, and (at large n) steal_count + load_imbalance if emitted.
#
# OUTPUTS
#   parallel_scaling.tsv          raw, one row per (n, threads, run)
#   parallel_scaling_summary.tsv  mean wall + END-TO-END speedup/efficiency
#   parallel_dp_speedup.tsv       DP-phase-only speedup/efficiency
#   amdahl_analysis.tsv           sequential-fraction estimate per graph size
#
# SUBMIT:  sbatch slurm/scripts/exp4_parallel_scaling.sh
# =============================================================================

#SBATCH --job-name=g2g_exp4
#SBATCH --output=slurm/logs/exp4-%j.out
#SBATCH --error=slurm/logs/exp4-%j.err
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=512GB
#SBATCH --constraint=cpu_amd_epyc_7763

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${G2G_REPO:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
BUILD="$REPO/build"
G2G="$BUILD/g2g"
GEN="$REPO/slurm/synth/gen_hprc_scale_gfa.py"
OUTDIR="$REPO/results/exp4"
GRAPHDIR="${G2G_GRAPHDIR:-$REPO/results/exp3/graphs}"   # reuse Exp3 graphs
RUNS="${G2G_RUNS:-5}"

# Graph sizes at which to run the thread sweep. We sweep at MULTIPLE sizes to
# show the bottleneck shifting: at 100K setup dominates; at 1M+ the DP dominates
# and the end-to-end curve approaches the DP curve. This is the core figure.
SCALE_SIZES="${G2G_SCALE_SIZES:-100000 1000000 10000000}"

THREAD_COUNTS="${G2G_THREADS:-1 2 4 8 16 32 64 128}"
NCPU="$(nproc)"

mkdir -p "$OUTDIR" "$GRAPHDIR" "$REPO/slurm/logs"
source "$SCRIPT_DIR/lib_measure.sh"

if [[ -n "${CONDA_BASE:-}" && -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
fi
conda activate "${G2G_ENV:-phase5}" 2>/dev/null || log "WARN: conda env not activated"

capture_env "$OUTDIR/env.txt"
log "exp4 start  sizes=[$SCALE_SIZES]  threads=[$THREAD_COUNTS]  runs=$RUNS"

# ── Build ────────────────────────────────────────────────────────────────────
if command -v cmake >/dev/null 2>&1; then
    cmake -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_CXX_COMPILER="$(which g++)" -DCMAKE_C_COMPILER="$(which gcc)" \
          -DCMAKE_CXX_FLAGS="-O3 -march=native" >/dev/null 2>&1 \
      && cmake --build "$BUILD" -j"$NCPU" --target g2g >/dev/null 2>&1 \
      && log "binary rebuilt" || log "WARN: rebuild failed; using existing binary"
fi
[[ -x "$G2G" ]] || { log "FATAL: no binary at $G2G"; exit 1; }

VERBOSE_FLAG=""
"$G2G" --help 2>&1 | grep -q -- "--verbose-timing" && VERBOSE_FLAG="--verbose-timing"
[[ -n "$VERBOSE_FLAG" ]] && log "per-phase timing ENABLED" || log "per-phase timing UNAVAILABLE (cols NA)"

# Does the binary print RunStats (tasks_executed / steal_count / imbalance)?
# Expected stderr (emit when --verbose-timing): [runstats] tasks=12345 steals=678 imbalance=1.23
parse_runstats() {
    local logfile="$1" line
    line=$(grep -m1 "\[runstats\]" "$logfile" 2>/dev/null || true)
    local tasks steals imb
    tasks=$( echo "$line" | grep -oP 'tasks=\K[0-9]+'         || echo NA)
    steals=$(echo "$line" | grep -oP 'steals=\K[0-9]+'        || echo NA)
    imb=$(   echo "$line" | grep -oP 'imbalance=\K[0-9.]+'    || echo NA)
    echo "${tasks:-NA} ${steals:-NA} ${imb:-NA}"
}
parse_phases() {
    local logfile="$1" line
    line=$(grep -m1 "\[timing\]" "$logfile" 2>/dev/null || true)
    [[ -z "$line" ]] && { echo "NA NA NA NA"; return; }
    local io idx dec dp
    io=$( echo "$line" | grep -oP 'io=\K[0-9.]+'        || echo NA)
    idx=$(echo "$line" | grep -oP 'index=\K[0-9.]+'     || echo NA)
    dec=$(echo "$line" | grep -oP 'decompose=\K[0-9.]+' || echo NA)
    dp=$( echo "$line" | grep -oP 'dp=\K[0-9.]+'        || echo NA)
    echo "${io:-NA} ${idx:-NA} ${dec:-NA} ${dp:-NA}"
}

RESULTS="$OUTDIR/parallel_scaling.tsv"
tsv_init "$RESULTS" "n_snarls\tthreads\trun\twall_sec\tpeak_rss_kb\tcpu_pct\texit\tphase_io\tphase_index\tphase_decompose\tphase_dp\ttasks_executed\tsteal_count\tload_imbalance"

# ── Sweep over graph sizes, then threads ─────────────────────────────────────
for N in $SCALE_SIZES; do
    GA="$GRAPHDIR/ga_${N}.gfa"; GB="$GRAPHDIR/gb_${N}.gfa"
    if [[ ! -s "$GA" || ! -s "$GB" ]]; then
        log "generating n=$N graphs..."
        python3 "$GEN" --mode size --n_snarls "$N" --out_dir "$GRAPHDIR" --seed 42
    fi

    # Warm the page cache once so the first thread count isn't penalised for I/O.
    log "warming cache for n=$N..."
    G2G_NUM_THREADS=1 "$G2G" "$GA" "$GB" --gaf /dev/null >/dev/null 2>&1 || true

    for T in $THREAD_COUNTS; do
        if (( T > NCPU )); then
            log "skip T=$T (> $NCPU CPUs)"; continue
        fi
        log "----- n=$N threads=$T -----"
        for RUN in $(seq 1 "$RUNS"); do
            g2glog=$(mktemp)
            # G2G_NUM_THREADS env hook makes this a no-recompile sweep.
            metrics=$(measure_full_capture "$g2glog" \
                env G2G_NUM_THREADS="$T" \
                "$G2G" "$GA" "$GB" --gaf /dev/null --threads "$T" $VERBOSE_FLAG)
            read -r wall peak user sys cpu ec <<< "$metrics"
            read -r p_io p_idx p_dec p_dp <<< "$(parse_phases "$g2glog")"
            read -r tasks steals imb <<< "$(parse_runstats "$g2glog")"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$N" "$T" "$RUN" "$wall" "$peak" "$cpu" "$ec" \
                "$p_io" "$p_idx" "$p_dec" "$p_dp" "$tasks" "$steals" "$imb" >> "$RESULTS"
            log "  T=$T run=$RUN wall=${wall}s dp=${p_dp}s tasks=$tasks steals=$steals imb=$imb"
            rm -f "$g2glog"
        done
    done
done

# ── Compute END-TO-END + DP-phase speedup, and Amdahl fraction ───────────────
log "computing speedup curves + Amdahl analysis..."
python3 - "$RESULTS" "$OUTDIR" <<'PYEOF'
import sys, csv, math
from collections import defaultdict
infile, outdir = sys.argv[1], sys.argv[2]

# rows[(n,threads)] = {"total":[...], "dp":[...], "tasks":set()}
rows = defaultdict(lambda: {"total": [], "dp": [], "tasks": set()})
with open(infile) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r.get("exit","0") not in ("0",""):
            continue
        key = (int(r["n_snarls"]), int(r["threads"]))
        try: rows[key]["total"].append(float(r["wall_sec"]))
        except: pass
        if r.get("phase_dp","NA") not in ("NA",""):
            try: rows[key]["dp"].append(float(r["phase_dp"]))
            except: pass
        if r.get("tasks_executed","NA") not in ("NA",""):
            rows[key]["tasks"].add(r["tasks_executed"])

def mean(v): return sum(v)/len(v) if v else 0.0
def stderr(v):
    if len(v)<2: return 0.0
    m=mean(v); return math.sqrt(sum((x-m)**2 for x in v)/(len(v)-1)/len(v))

sizes = sorted({n for (n,_) in rows})

# END-TO-END summary
with open(f"{outdir}/parallel_scaling_summary.tsv","w") as f:
    f.write("n_snarls\tthreads\tmean_wall_sec\tstderr_wall_sec\tspeedup\tefficiency\ttasks_const\tn_runs\n")
    for n in sizes:
        t1=mean(rows[(n,1)]["total"]) if (n,1) in rows else 0.0
        for (nn,t) in sorted(k for k in rows if k[0]==n):
            m=mean(rows[(nn,t)]["total"]); se=stderr(rows[(nn,t)]["total"])
            sp = t1/m if m>0 and t1>0 else 0.0
            eff= sp/t if t>0 else 0.0
            tasks = rows[(nn,t)]["tasks"]
            const = "yes" if len(tasks)<=1 else f"NO({len(tasks)} distinct)"
            f.write(f"{nn}\t{t}\t{m:.4f}\t{se:.4f}\t{sp:.3f}\t{eff:.3f}\t{const}\t{len(rows[(nn,t)]['total'])}\n")

# DP-PHASE-ONLY summary
with open(f"{outdir}/parallel_dp_speedup.tsv","w") as f:
    f.write("n_snarls\tthreads\tmean_dp_sec\tdp_speedup\tdp_efficiency\tn_runs\n")
    for n in sizes:
        dp1=mean(rows[(n,1)]["dp"]) if (n,1) in rows else 0.0
        for (nn,t) in sorted(k for k in rows if k[0]==n):
            mdp=mean(rows[(nn,t)]["dp"])
            sp = dp1/mdp if mdp>0 and dp1>0 else 0.0
            eff= sp/t if t>0 else 0.0
            f.write(f"{nn}\t{t}\t{mdp:.4f}\t{sp:.3f}\t{eff:.3f}\t{len(rows[(nn,t)]['dp'])}\n")

# AMDAHL: estimate sequential fraction s from 1-thread phase split.
# s ~ (io+index+decompose) / total ; parallelisable p = dp/total.
# Max theoretical speedup = 1/(s + (1-s)/P).
with open(f"{outdir}/amdahl_analysis.tsv","w") as f:
    f.write("n_snarls\tmean_total_1t\tmean_dp_1t\tseq_fraction\tpar_fraction\tmax_speedup_128\n")
    for n in sizes:
        if (n,1) not in rows: continue
        tot=mean(rows[(n,1)]["total"]); dp=mean(rows[(n,1)]["dp"])
        if tot<=0: continue
        par = dp/tot if dp>0 else 0.0
        s = max(0.0, 1.0-par)
        max_sp = 1.0/(s + (1-s)/128) if (s+(1-s)/128)>0 else 0.0
        f.write(f"{n}\t{tot:.4f}\t{dp:.4f}\t{s:.4f}\t{par:.4f}\t{max_sp:.2f}\n")

print("[exp4] wrote parallel_scaling_summary.tsv, parallel_dp_speedup.tsv, amdahl_analysis.tsv")
PYEOF

log "=== Exp4 DONE ==="
log "  raw:           $RESULTS"
log "  end-to-end:    $OUTDIR/parallel_scaling_summary.tsv"
log "  dp-only:       $OUTDIR/parallel_dp_speedup.tsv"
log "  amdahl:        $OUTDIR/amdahl_analysis.tsv"
