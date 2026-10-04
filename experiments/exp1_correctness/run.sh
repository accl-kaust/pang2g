#!/usr/bin/env bash
# ============================================================
# Experiment 1 — Synthetic Ground-Truth Correctness
#
# Runs the Exp1 test suite (planted variants
# at known positions; measures recall and precision).
# Then runs the pairing-error sweep (Exp1b) across graph sizes to
# produce the precision/recall data behind fig_pairing_error.
#
# Usage:
#   conda activate g2g_eval
#   bash experiments/exp1_correctness/run.sh
#
# Output: results/phase5/exp1/
# ============================================================
set -euo pipefail

G2G_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RESULTS="$G2G_ROOT/results/phase5/exp1"
BUILD="$G2G_ROOT/build/release"
BINARY="$BUILD/tests/test_phase5_exp1"
TIMESTAMP=$(date +%Y%m%dT%H%M%S)

mkdir -p "$RESULTS"

echo "[exp1] G2G root:   $G2G_ROOT"
echo "[exp1] Results:    $RESULTS"
echo "[exp1] Timestamp:  $TIMESTAMP"

# ─────────────────────────────────────────────────────────────
# 0. Build if needed
# ─────────────────────────────────────────────────────────────
if [[ ! -f "$BINARY" ]]; then
    echo "[exp1] Binary not found — building..."
    cmake -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -G Ninja "$G2G_ROOT"
    cmake --build "$BUILD" -j"$(nproc)"
fi

# ─────────────────────────────────────────────────────────────
# 1. Run the ctest suite (43 + standalone tests)
# ─────────────────────────────────────────────────────────────
echo "[exp1] Running ctest suite..."
START_NS=$(date +%s%N)

cd "$BUILD"
ctest --output-on-failure \
      -R "phase5_exp1" \
      --output-junit "$RESULTS/junit_exp1_${TIMESTAMP}.xml" \
      2>&1 | tee "$RESULTS/ctest_log_${TIMESTAMP}.txt"

END_NS=$(date +%s%N)
ELAPSED_MS=$(( (END_NS - START_NS) / 1000000 ))
echo "[exp1] ctest suite finished in ${ELAPSED_MS} ms"

# ─────────────────────────────────────────────────────────────
# 2. Pairing-error sweep (Exp1b): planted graphs of 5 to 1e5 snarls, short
#    (colliding) boundary keys and private snarls, boundary-key vs.
#    anchor-interval pairing.  Scores BOTH deltas by allele sequence.
#    (The earlier inline sweep read the wrong VCF column and is removed.)
# ─────────────────────────────────────────────────────────────
echo "[exp1] Running pairing-error sweep..."
SWEEP_CSV="$RESULTS/exp1b_pairing_error_${TIMESTAMP}.csv"
python3 "$G2G_ROOT/experiments/exp1_correctness/pairing_error_sweep.py" \
    --g2g "$BUILD/g2g" \
    --pairing both --anchor-lens 8,2,1 --private-fracs 0,0.1 --linker-every 4 \
    --out "$SWEEP_CSV"

echo "[exp1] Sweep results: $SWEEP_CSV"
echo "[exp1] Figure: python3 make_figures.py --only fig_pairing_error --pairing-csv $SWEEP_CSV"
echo "[exp1] Done."
