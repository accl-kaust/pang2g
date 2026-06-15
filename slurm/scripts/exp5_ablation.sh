#!/bin/bash --login
# =============================================================================
# Experiment 5 — Ablation Study (Q3: Does the system design matter?)
# =============================================================================
# Roadmap §2.7.3 / Phase 5, Experiment 5
#
# Five ablations measured on a fixed graph size (10K snarls):
#
#   A1. Cache-oblivious (vEB) layout vs. BFS layout
#       → metric: L3 miss rate via perf stat; throughput
#       → implementation: g2g built with G2G_DISABLE_VEB=1 compile flag
#         (add this flag to decompose.cpp to skip apply_veb_order)
#
#   A2. Work-stealing (steal-largest) vs. static round-robin partition
#       → metric: wall time; load imbalance (max_task_time / mean_task_time)
#       → implementation: G2G_STEAL_RANDOM=1 env flag (steal random instead
#         of priority queue); static partition = --sequential with
#         omp-style manual chunking (approximated by single thread)
#
#   A3. Steal-largest vs. steal-random
#       → metric: wall time on heavy-tailed graph (Zipf alpha=1.0 → very skewed)
#
#   A4. MinimizerIndex vs. naive O(|VA||VB|) anchor detection
#       → metric: anchor detection wall time (logged by g2g stderr)
#       → method: run on graph where |VA|x|VB| > 10^6 triggers the index path;
#         compare vs patched binary that always uses naive path
#
#   A5. Hungarian (LAPJV) vs. brute-force matching
#       → metric: matching time at k=5,10,20,50
#       → implementation: synthetic micro-benchmark (matching_bench binary)
#
# NOTE: A1 (perf stat L3 miss rate) requires perf to be available.
#       On the SANDS testbed this needs:
#         sysctl kernel.perf_event_paranoid=2  (ask admin or use sudo)
#       The script checks for perf and gracefully skips A1 perf counters
#       if unavailable, falling back to wall-time only.
#
# Outputs
#   ablation_results.tsv    — all ablation measurements
#
# Usage
#   sbatch slurm/scripts/exp5_ablation.sh
# =============================================================================

#SBATCH --job-name=g2g_exp5_ablation
#SBATCH --output=slurm/logs/exp5-%j.out
#SBATCH --error=slurm/logs/exp5-%j.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=128GB
#SBATCH --constraint=cpu_amd_epyc_7763

set -euo pipefail

REPO="$HOME/Research/G2G/G2G_Project_Summer_2026"
BUILD_BASE="$REPO/build"
SYNTH_GEN="$REPO/slurm/synth/gen_scaling_gfa.py"
OUTDIR="$REPO/results/exp5"
GRAPHDIR="$REPO/results/exp3/graphs"
LOGDIR="$REPO/slurm/logs"
RUNS=5

mkdir -p "$OUTDIR" "$GRAPHDIR" "$LOGDIR"
source "$HOME/miniconda3/etc/profile.d/conda.sh" && conda activate phase5

RESULTS="$OUTDIR/ablation_results.tsv"
echo -e "ablation\tvariant\tparam\trun\twall_sec\tl3_miss_rate" > "$RESULTS"

# ── Wall-time measurement helper ─────────────────────────────────────────────
wall_sec() {
    local tmp; tmp=$(mktemp)
    /usr/bin/time -v "$@" 2>"$tmp"
    grep "Elapsed (wall clock)" "$tmp" | awk '{print $NF}' | \
        awk -F: '{if(NF==3) print $1*3600+$2*60+$3; else print $1*60+$2}'
    rm -f "$tmp"
}

# ── perf stat helper (L3 cache miss rate) ────────────────────────────────────
HAVE_PERF=0
if perf stat true 2>/dev/null; then HAVE_PERF=1; fi

perf_l3_miss_rate() {
    # Returns L3-miss / L3-reference ratio (0.0–1.0), or "N/A" if no perf.
    if [[ $HAVE_PERF -eq 0 ]]; then echo "N/A"; return; fi
    local tmp; tmp=$(mktemp)
    perf stat -e cache-references,cache-misses "$@" 2>"$tmp" >/dev/null
    local refs misses
    refs=$(  grep "cache-references" "$tmp" | awk '{gsub(",","",$1); print $1}')
    misses=$(grep "cache-misses"     "$tmp" | awk '{gsub(",","",$1); print $1}')
    rm -f "$tmp"
    if [[ -z "$refs" || "$refs" -eq 0 ]]; then echo "N/A"; else
        python3 -c "print(f'{$misses/$refs:.4f}')"
    fi
}

# ── Build helper: compile g2g with optional extra flags ──────────────────────
build_variant() {
    local suffix="$1"; shift   # e.g. "veb" or "noveb"
    local extra_flags="${1:-}"; shift || true
    local bdir="$BUILD_BASE/$suffix"
    cmake -B "$bdir" -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_CXX_COMPILER="$(which g++)" \
          -DCMAKE_C_COMPILER="$(which gcc)" \
          -DCMAKE_CXX_FLAGS="-O3 -march=native ${extra_flags}" \
          "$REPO" > /dev/null 2>&1
    cmake --build "$bdir" -j"$(nproc)" --target g2g > /dev/null 2>&1
    echo "$bdir/g2g"
}

# ── Generate graphs ───────────────────────────────────────────────────────────
N_STD=10000          # standard size for A1–A3
N_ANCHOR=200000      # size to trigger MinimizerIndex (|VA|x|VB| > 1e6)
N_HEAVY=10000        # heavy-tailed graph for A3 (Zipf alpha=1.0)

for N in $N_STD $N_ANCHOR; do
    GFA="$GRAPHDIR/ga_${N}.gfa"
    [[ -f "$GFA" ]] || python3 "$SYNTH_GEN" --n_snarls "$N" \
                               --out_dir "$GRAPHDIR" --seed 42
done

# Heavy-tailed graph: re-generate with seed=999 to get different distribution
GFA_HEAVY_A="$GRAPHDIR/ga_heavy_${N_HEAVY}.gfa"
GFA_HEAVY_B="$GRAPHDIR/gb_heavy_${N_HEAVY}.gfa"
if [[ ! -f "$GFA_HEAVY_A" ]]; then
    python3 "$SYNTH_GEN" --n_snarls "$N_HEAVY" \
            --out_dir "$GRAPHDIR" --seed 999 \
            && mv "$GRAPHDIR/ga_${N_HEAVY}.gfa" "$GFA_HEAVY_A" \
            && mv "$GRAPHDIR/gb_${N_HEAVY}.gfa" "$GFA_HEAVY_B"
fi

GFA_A_STD="$GRAPHDIR/ga_${N_STD}.gfa"
GFA_B_STD="$GRAPHDIR/gb_${N_STD}.gfa"

# =============================================================================
# A1: vEB layout vs. BFS layout
# =============================================================================
echo "[exp5] === A1: vEB vs BFS layout ==="

# Standard build (vEB enabled)
G2G_VEB=$(build_variant "veb")
# Build with vEB disabled (DG2G_DISABLE_VEB preprocessor flag)
G2G_BFS=$(build_variant "bfs" "-DG2G_DISABLE_VEB=1")

for RUN in $(seq 1 "$RUNS"); do
    for VARIANT in veb bfs; do
        [[ "$VARIANT" == "veb" ]] && BIN="$G2G_VEB" || BIN="$G2G_BFS"
        # Sequential run to isolate layout effect from parallelism
        WALL=$(wall_sec "$BIN" "$GFA_A_STD" "$GFA_B_STD" \
               --sequential --gaf /dev/null)
        L3=$(perf_l3_miss_rate "$BIN" "$GFA_A_STD" "$GFA_B_STD" \
             --sequential --gaf /dev/null)
        echo -e "A1\t${VARIANT}\t${N_STD}\t${RUN}\t${WALL}\t${L3}" >> "$RESULTS"
        echo "[exp5]   A1 variant=$VARIANT run=$RUN wall=${WALL}s l3_miss=${L3}"
    done
done

# =============================================================================
# A2 & A3: Work-stealing steal-largest vs. steal-random
# =============================================================================
# steal-random approximation: build with G2G_STEAL_RANDOM=1 which replaces
# the priority queue pop with a random selection from the ready queue.
echo "[exp5] === A2/A3: steal-largest vs. steal-random ==="

G2G_SL=$(build_variant "steal_largest")                        # default
G2G_SR=$(build_variant "steal_random" "-DG2G_STEAL_RANDOM=1")  # random steal

for GRAPH_TAG in std heavy; do
    if [[ "$GRAPH_TAG" == "std" ]]; then
        FA="$GFA_A_STD"; FB="$GFA_B_STD"; PARAM="$N_STD"
    else
        FA="$GFA_HEAVY_A"; FB="$GFA_HEAVY_B"; PARAM="${N_HEAVY}_heavy"
    fi

    for RUN in $(seq 1 "$RUNS"); do
        for VARIANT in steal_largest steal_random; do
            [[ "$VARIANT" == "steal_largest" ]] && BIN="$G2G_SL" || BIN="$G2G_SR"
            WALL=$(wall_sec "$BIN" "$FA" "$FB" \
                   --threads "$(nproc)" --gaf /dev/null)
            echo -e "A2_A3\t${VARIANT}\t${PARAM}\t${RUN}\t${WALL}\tN/A" >> "$RESULTS"
            echo "[exp5]   A2/A3 variant=$VARIANT graph=$GRAPH_TAG run=$RUN wall=${WALL}s"
        done
    done
done

# =============================================================================
# A4: MinimizerIndex vs. naive anchor detection
# =============================================================================
# Use a large graph where the index path is triggered automatically.
# Build a "naive-only" variant with G2G_DISABLE_MINIMIZER_INDEX=1.
echo "[exp5] === A4: MinimizerIndex vs. naive ==="

GFA_A_ANC="$GRAPHDIR/ga_${N_ANCHOR}.gfa"
GFA_B_ANC="$GRAPHDIR/gb_${N_ANCHOR}.gfa"

G2G_IDX=$(build_variant "minimizer_idx")
G2G_NAIVE=$(build_variant "naive_anchor" "-DG2G_DISABLE_MINIMIZER_INDEX=1")

for RUN in $(seq 1 "$RUNS"); do
    for VARIANT in minimizer_idx naive_anchor; do
        [[ "$VARIANT" == "minimizer_idx" ]] && BIN="$G2G_IDX" || BIN="$G2G_NAIVE"
        WALL=$(wall_sec "$BIN" "$GFA_A_ANC" "$GFA_B_ANC" \
               --sequential --gaf /dev/null)
        echo -e "A4\t${VARIANT}\t${N_ANCHOR}\t${RUN}\t${WALL}\tN/A" >> "$RESULTS"
        echo "[exp5]   A4 variant=$VARIANT run=$RUN wall=${WALL}s"
    done
done

# =============================================================================
# A5: Hungarian (LAPJV) vs. brute-force — matching micro-benchmark
# =============================================================================
# Write and build a tiny standalone benchmark that exercises matching alone.
echo "[exp5] === A5: Hungarian vs. brute-force matching ==="

BENCH_SRC="$REPO/slurm/synth/matching_bench.cpp"
BENCH_BIN="$OUTDIR/matching_bench"

cat > "$BENCH_SRC" <<'CPPSRC'
/**
 * matching_bench.cpp — micro-benchmark for allele_match dispatch.
 * Usage: matching_bench <k> <n_trials>
 */
#include <g2g/matching.hpp>
#include <g2g/sequence.hpp>
#include <chrono>
#include <cstdlib>
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
    if (argc < 3) { std::cerr << "Usage: matching_bench <k> <n_trials>\n"; return 1; }
    int k        = std::atoi(argv[1]);
    int n_trials = std::atoi(argv[2]);

    std::mt19937 rng(42);
    AlignerParams params;

    // Build k x k similarity matrix
    std::vector<std::string> sa(k), sb(k);
    for (auto& s : sa) s = rand_dna(50, rng);
    for (auto& s : sb) s = rand_dna(50, rng);

    std::vector<std::vector<float>> sim(k, std::vector<float>(k));
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j)
            sim[i][j] = seq_similarity(sa[i], sb[j]);

    // Warm up
    for (int t = 0; t < 10; ++t) allele_match(sim, params);

    auto t0 = Clock::now();
    for (int t = 0; t < n_trials; ++t) allele_match(sim, params);
    auto t1 = Clock::now();
    double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
    std::cout << k << "\t" << (us / n_trials) << "\n";   // k, microseconds/trial
    return 0;
}
CPPSRC

# Compile micro-benchmark
g++ -std=c++20 -O3 -march=native \
    -I"$REPO/include" \
    "$BENCH_SRC" \
    "$BUILD_BASE/steal_largest/libg2g_lib.a" \
    -pthread -fopenmp \
    -o "$BENCH_BIN"

BENCH_RESULTS="$OUTDIR/matching_bench.tsv"
echo -e "k\tmicros_per_trial\trun" > "$BENCH_RESULTS"

for K in 2 3 5 7 10 20 50; do
    for RUN in $(seq 1 "$RUNS"); do
        LINE=$("$BENCH_BIN" "$K" 1000)
        echo -e "${LINE}\t${RUN}" >> "$BENCH_RESULTS"
        echo "[exp5]   A5 k=$K run=$RUN: ${LINE} µs/trial"
    done
done

echo "[exp5] Results:         $RESULTS"
echo "[exp5] Matching bench:  $BENCH_RESULTS"
echo "[exp5] DONE"
