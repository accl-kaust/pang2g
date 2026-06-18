#!/usr/bin/env bash
# ============================================================
# Experiment 1 — Synthetic Ground-Truth Correctness
#
# Runs the Phase 5 Exp1 test suite (43 tests: planted variants
# at known positions; measures recall and precision).
# Also runs a sweep across (n_snarls, k_shared, k_only) to
# produce a recall/precision table for the paper (Figure 1).
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
# 2. Extended parameter sweep using the synthetic GFA generator
#    Varies: n_snarls, k_alleles, shared_fraction
#    Uses the g2g binary directly on GFA files.
# ─────────────────────────────────────────────────────────────
echo "[exp1] Running parameter sweep..."

G2G_BIN="$BUILD/g2g"
GEN_SCRIPT="$G2G_ROOT/experiments/exp3_throughput/gen_synthetic_gfa.py"
SWEEP_CSV="$RESULTS/exp1_recall_precision_sweep_${TIMESTAMP}.csv"

echo "n_snarls,k_shared,k_only,seed,recall_a,precision_a,recall_b,precision_b,f1_a,f1_b,wall_ms,peak_kb" \
    > "$SWEEP_CSV"

N_SNARLS_LIST=(5 10 20 50 100 200)
K_SHARED_LIST=(1 2 3)
K_ONLY_LIST=(1 2)
SEEDS=(42 123 999)

for n in "${N_SNARLS_LIST[@]}"; do
for ks in "${K_SHARED_LIST[@]}"; do
for ko in "${K_ONLY_LIST[@]}"; do
for seed in "${SEEDS[@]}"; do
    TMPDIR_RUN=$(mktemp -d "$RESULTS/run_XXXXXX")
    trap "rm -rf $TMPDIR_RUN" EXIT

    # Generate GFA
    python3 "$GEN_SCRIPT" \
        --mode planted \
        --n-snarls "$n" \
        --k-shared "$ks" \
        --k-only "$ko" \
        --seed "$seed" \
        --out-a "$TMPDIR_RUN/ga.gfa" \
        --out-b "$TMPDIR_RUN/gb.gfa" \
        --truth-a "$TMPDIR_RUN/truth_a.txt" \
        --truth-b "$TMPDIR_RUN/truth_b.txt"

    # Run g2g with timing
    START_NS=$(date +%s%N)
    /usr/bin/time -v -o "$TMPDIR_RUN/time.txt" \
        "$G2G_BIN" "$TMPDIR_RUN/ga.gfa" "$TMPDIR_RUN/gb.gfa" \
        --gaf "$TMPDIR_RUN/out.gaf" \
        --vcf "$TMPDIR_RUN/out.vcf" \
        2>"$TMPDIR_RUN/g2g_stderr.txt"
    END_NS=$(date +%s%N)
    WALL_MS=$(( (END_NS - START_NS) / 1000000 ))

    # Peak memory (KB)
    PEAK_KB=$(grep "Maximum resident" "$TMPDIR_RUN/time.txt" 2>/dev/null | \
              grep -o '[0-9]*' | tail -1 || echo 0)

    # Compute recall/precision from VCF and ground truth
    python3 - <<PYEOF
import sys
def load_truth(f):
    try:
        with open(f) as fh:
            return set(l.strip() for l in fh if l.strip())
    except FileNotFoundError:
        return set()

def load_reported(vcf):
    seqs = set()
    try:
        with open(vcf) as fh:
            for line in fh:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) > 7:
                    seqs.add(parts[7].split(';')[0].replace('SEQ=',''))
    except FileNotFoundError:
        pass
    return seqs

truth_a = load_truth("$TMPDIR_RUN/truth_a.txt")
truth_b = load_truth("$TMPDIR_RUN/truth_b.txt")
rep_a   = load_reported("$TMPDIR_RUN/out.vcf")

tp_a = len(truth_a & rep_a)
fp_a = len(rep_a - truth_a)
fn_a = len(truth_a - rep_a)
prec_a = tp_a/(tp_a+fp_a) if tp_a+fp_a>0 else 1.0
rec_a  = tp_a/(tp_a+fn_a) if tp_a+fn_a>0 else 1.0
f1_a   = 2*prec_a*rec_a/(prec_a+rec_a) if prec_a+rec_a>0 else 0.0

print(f"$n,$ks,$ko,$seed,{rec_a:.4f},{prec_a:.4f},0.0,0.0,{f1_a:.4f},0.0,$WALL_MS,$PEAK_KB")
PYEOF

done; done; done; done >> "$SWEEP_CSV"

echo "[exp1] Sweep results: $SWEEP_CSV"
echo "[exp1] Done."
