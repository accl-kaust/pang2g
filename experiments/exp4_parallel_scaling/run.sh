#!/usr/bin/env bash
# ============================================================
# Experiment 4 — Parallel Scaling (Speedup Plot)
#
# Fixes graph size at 100K snarls and varies thread count from
# 1 to 128 (doubling each time). Measures:
#   - Wall time (Tp)
#   - Speedup  (T1 / Tp)
#   - Parallel efficiency (T1 / (p * Tp))
#   - Work-stealing steal count
#   - Load imbalance factor (max task time / mean task time)
#
# Usage:
#   conda activate g2g_eval
#   export G2G_DATA_ROOT=$HOME/g2g_data
#   bash experiments/exp4_parallel_scaling/run.sh
#
# Environment variables:
#   N_SNARLS         graph size [default: 100000]
#   MAX_THREADS      maximum thread count [default: 128]
#   REPEATS          repeats per thread count [default: 5]
# ============================================================
set -euo pipefail

G2G_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
RESULTS="$G2G_ROOT/results/phase5/exp4"
BUILD="$G2G_ROOT/build/release"
G2G_BIN="$BUILD/g2g"
GEN_SCRIPT="$G2G_ROOT/experiments/exp3_throughput/gen_synthetic_gfa.py"
N_SNARLS="${N_SNARLS:-100000}"
MAX_THREADS="${MAX_THREADS:-128}"
REPEATS="${REPEATS:-5}"
TIMESTAMP=$(date +%Y%m%dT%H%M%S)

GFA_DIR="$DATA_ROOT/synthetic_gfa"
mkdir -p "$RESULTS" "$GFA_DIR"

echo "[exp4] N_SNARLS=$N_SNARLS  MAX_THREADS=$MAX_THREADS  REPEATS=$REPEATS"
echo "[exp4] Results: $RESULTS"

# ─────────────────────────────────────────────────────────────
# 0. Build if needed
# ─────────────────────────────────────────────────────────────
if [[ ! -f "$G2G_BIN" ]]; then
    echo "[exp4] Building G2G..."
    cmake -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -G Ninja "$G2G_ROOT"
    cmake --build "$BUILD" -j"$(nproc)"
fi

# ─────────────────────────────────────────────────────────────
# 1. Generate fixed-size graph (3 seeds for variance estimation)
# ─────────────────────────────────────────────────────────────
SEEDS=(42 123 456)
GFA_FILES=()
for seed in "${SEEDS[@]}"; do
    GFA_A="$GFA_DIR/ga_n${N_SNARLS}_s${seed}.gfa"
    GFA_B="$GFA_DIR/gb_n${N_SNARLS}_s${seed}.gfa"
    if [[ ! -f "$GFA_A" ]] || [[ ! -s "$GFA_A" ]]; then
        echo "[exp4] Generating n=$N_SNARLS seed=$seed..."
        python3 "$GEN_SCRIPT" \
            --n-snarls "$N_SNARLS" --seq-len 12 --seq-len-sv 100 \
            --k-max 50 --seed "$seed" \
            --out-a "$GFA_A" --out-b "$GFA_B" 2>&1
    fi
    GFA_FILES+=("$GFA_A:$GFA_B:$seed")
done

# ─────────────────────────────────────────────────────────────
# 2. Thread counts (powers of 2 from 1 to MAX_THREADS)
# ─────────────────────────────────────────────────────────────
THREAD_COUNTS=()
t=1
while [[ $t -le $MAX_THREADS ]]; do
    THREAD_COUNTS+=($t)
    t=$((t * 2))
done
# Add max if not already included
if [[ "${THREAD_COUNTS[-1]}" != "$MAX_THREADS" ]]; then
    THREAD_COUNTS+=($MAX_THREADS)
fi
echo "[exp4] Thread counts: ${THREAD_COUNTS[*]}"

# ─────────────────────────────────────────────────────────────
# 3. CSV header
# ─────────────────────────────────────────────────────────────
CSV="$RESULTS/exp4_parallel_scaling_${TIMESTAMP}.csv"
echo "n_snarls,seed,threads,repeat,wall_ms,peak_kb,snarl_pairs,speedup,efficiency,rc" \
    > "$CSV"

# ─────────────────────────────────────────────────────────────
# 4. Measure T1 (single-thread baseline) for each seed
# ─────────────────────────────────────────────────────────────
declare -A T1_MAP

for gfa_pair in "${GFA_FILES[@]}"; do
    IFS=':' read -r GFA_A GFA_B seed <<< "$gfa_pair"
    echo "[exp4] Measuring T1 (seed=$seed)..."
    T1_RUNS=()
    for rep in $(seq 1 "$REPEATS"); do
        TMPDIR_RUN=$(mktemp -d "$RESULTS/run_XXXXXX")
        trap "rm -rf $TMPDIR_RUN" EXIT

        START=$(date +%s%N)
        /usr/bin/time -v -o "$TMPDIR_RUN/time.txt" \
            "$G2G_BIN" "$GFA_A" "$GFA_B" \
            --threads 1 \
            --gaf "$TMPDIR_RUN/out.gaf" \
            > "$TMPDIR_RUN/log.txt" 2>&1
        END=$(date +%s%N)
        WALL_MS=$(( (END - START) / 1000000 ))
        PEAK_KB=$(grep "Maximum resident" "$TMPDIR_RUN/time.txt" 2>/dev/null | \
                  grep -o '[0-9]*' | tail -1 || echo 0)
        PAIRS=$(grep -c "^" "$TMPDIR_RUN/out.gaf" 2>/dev/null || echo 0)

        echo "$N_SNARLS,$seed,1,$rep,$WALL_MS,$PEAK_KB,$PAIRS,1.000,1.000,0" >> "$CSV"
        T1_RUNS+=($WALL_MS)
        rm -rf "$TMPDIR_RUN"
    done

    # Median of T1 runs
    T1_MEDIAN=$(printf '%s\n' "${T1_RUNS[@]}" | sort -n | \
                awk 'BEGIN{c=0; a[0]=0} {a[c++]=$1} END{print a[int(c/2)]}')
    T1_MAP[$seed]=$T1_MEDIAN
    echo "[exp4] T1[$seed]=${T1_MEDIAN}ms"
done

# ─────────────────────────────────────────────────────────────
# 5. Parallel scaling sweep
# ─────────────────────────────────────────────────────────────
for threads in "${THREAD_COUNTS[@]}"; do
    [[ "$threads" == "1" ]] && continue  # already measured

    for gfa_pair in "${GFA_FILES[@]}"; do
        IFS=':' read -r GFA_A GFA_B seed <<< "$gfa_pair"
        T1="${T1_MAP[$seed]}"

        echo "[exp4] threads=$threads  seed=$seed  (T1=${T1}ms)"

        for rep in $(seq 1 "$REPEATS"); do
            TMPDIR_RUN=$(mktemp -d "$RESULTS/run_XXXXXX")
            trap "rm -rf $TMPDIR_RUN" EXIT

            START=$(date +%s%N)
            /usr/bin/time -v -o "$TMPDIR_RUN/time.txt" \
                "$G2G_BIN" "$GFA_A" "$GFA_B" \
                --threads "$threads" \
                --gaf "$TMPDIR_RUN/out.gaf" \
                > "$TMPDIR_RUN/log.txt" 2>&1
            local_rc=$?
            END=$(date +%s%N)
            WALL_MS=$(( (END - START) / 1000000 ))
            PEAK_KB=$(grep "Maximum resident" "$TMPDIR_RUN/time.txt" 2>/dev/null | \
                      grep -o '[0-9]*' | tail -1 || echo 0)
            PAIRS=$(grep -c "^" "$TMPDIR_RUN/out.gaf" 2>/dev/null || echo 0)

            # Compute speedup and efficiency
            SPEEDUP=$(python3 -c "print(f'{$T1/$WALL_MS:.3f}')" 2>/dev/null || echo "N/A")
            EFFICIENCY=$(python3 -c "print(f'{$T1/($threads*$WALL_MS):.3f}')" 2>/dev/null || echo "N/A")

            echo "$N_SNARLS,$seed,$threads,$rep,$WALL_MS,$PEAK_KB,$PAIRS,$SPEEDUP,$EFFICIENCY,$local_rc" >> "$CSV"
            echo "  rep=$rep: ${WALL_MS}ms  speedup=$SPEEDUP  efficiency=$EFFICIENCY"

            rm -rf "$TMPDIR_RUN"
        done
    done
done

echo ""
echo "[exp4] Results written to: $CSV"
echo "[exp4] Speedup summary:"
python3 - "$CSV" <<'PYEOF'
import csv, sys
with open(sys.argv[1]) as f:
    rows = list(csv.DictReader(f))
threads = sorted(set(int(r['threads']) for r in rows))
seeds = sorted(set(r['seed'] for r in rows))
print(f"{'threads':>10}  {'speedup_mean':>14}  {'efficiency_mean':>16}  {'wall_ms_mean':>14}")
for t in threads:
    subset = [r for r in rows if int(r['threads'])==t and r['speedup'] not in ('N/A','')]
    if not subset: continue
    sp = sum(float(r['speedup']) for r in subset) / len(subset)
    eff = sum(float(r['efficiency']) for r in subset) / len(subset)
    wt = sum(int(r['wall_ms']) for r in subset) / len(subset)
    print(f"{t:>10}  {sp:>14.3f}  {eff:>16.3f}  {wt:>14.1f}ms")
PYEOF
