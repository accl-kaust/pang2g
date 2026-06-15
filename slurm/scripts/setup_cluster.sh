#!/bin/bash
# =============================================================================
# setup_cluster.sh — One-time setup on mcmgt01 before running Phase 5 jobs.
#
# Run this interactively (NOT via sbatch):
#   bash ~/Research/G2G/G2G_Project_Summer_2026/slurm/scripts/setup_cluster.sh
#
# Assumptions:
#   - Repo is already cloned at ~/Research/G2G/G2G_Project_Summer_2026
#   - Miniconda3 is installed at ~/miniconda3
# =============================================================================

set -euo pipefail

REPO="$HOME/Research/G2G/G2G_Project_Summer_2026"
CONDA_BASE="$HOME/miniconda3"
ENV_NAME="phase5"

echo "======================================================="
echo "  G2G Phase 5 — Cluster Setup"
echo "  Node: $(hostname)"
echo "  Date: $(date)"
echo "  Repo: $REPO"
echo "======================================================="

# ── 1. Initialise conda from miniconda3 ──────────────────────────────────────
echo "[setup] Initialising conda from $CONDA_BASE..."
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate base
echo "[setup] conda: $(conda --version)"

# ── 2. Create or update the phase5 environment ───────────────────────────────
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "[setup] Environment '${ENV_NAME}' exists — updating packages..."
    conda activate "$ENV_NAME"
    conda install -y -c conda-forge -c bioconda cmake cxx-compiler minimap2
else
    echo "[setup] Creating environment '${ENV_NAME}'..."
    conda create -y -n "$ENV_NAME" \
        -c conda-forge -c bioconda \
        cmake cxx-compiler minimap2 python=3.11
    conda activate "$ENV_NAME"
fi

echo "[setup] Compiler:  $(g++ --version | head -1)"
echo "[setup] CMake:     $(cmake --version | head -1)"
echo "[setup] minimap2:  $(minimap2 --version)"
echo "[setup] OpenMP:    $(g++ -fopenmp -dM -E - < /dev/null 2>/dev/null | grep _OPENMP || echo 'not found')"

# ── 3. Build G2G (Release) ───────────────────────────────────────────────────
echo ""
echo "[setup] Building G2G (Release)..."
cmake -B "$REPO/build" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER="$(which g++)" \
      -DCMAKE_C_COMPILER="$(which gcc)" \
      -DCMAKE_CXX_FLAGS="-O3 -march=native"
cmake --build "$REPO/build" -j"$(nproc)"

# ── 4. Run full test suite ────────────────────────────────────────────────────
echo ""
echo "[setup] Running test suite..."
cd "$REPO/build"
ctest --output-on-failure -j"$(nproc)" 2>&1 | tail -12
cd "$REPO"

# ── 5. Create output and log directories ─────────────────────────────────────
mkdir -p "$REPO/results"/{exp3,exp4,exp5} "$REPO/slurm/logs"
echo "[setup] Output dirs: $REPO/results/{exp3,exp4,exp5}"

# ── 6. Check perf availability (needed for Exp5 A1 L3 miss-rate) ─────────────
echo ""
echo "[setup] Checking perf availability..."
if perf stat true 2>/dev/null; then
    echo "[setup] perf OK — L3 miss-rate will be collected in Exp5 A1"
else
    PARANOID=$(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo "unknown")
    echo "[setup] WARNING: perf unavailable (perf_event_paranoid=$PARANOID)"
    echo "[setup]   Exp5 A1 will report wall-time only (no L3 miss rate)."
    echo "[setup]   To fix: ask admin → sudo sysctl kernel.perf_event_paranoid=2"
fi

# ── 7. Storage check ─────────────────────────────────────────────────────────
echo ""
echo "[setup] Storage:"
df -h "$HOME" | tail -1

# ── 8. Smoke test ────────────────────────────────────────────────────────────
echo ""
echo "[setup] Smoke test..."
python3 "$REPO/slurm/synth/gen_scaling_gfa.py" \
    --n_snarls 100 --out_dir /tmp/g2g_smoke --seed 1 2>/dev/null
"$REPO/build/g2g" /tmp/g2g_smoke/ga_100.gfa /tmp/g2g_smoke/gb_100.gfa \
    --gaf /dev/null --threads 4 2>&1 | tail -3
rm -rf /tmp/g2g_smoke
echo "[setup] Smoke test PASSED"

echo ""
echo "======================================================="
echo "  Setup complete. Submit jobs from $REPO:"
echo ""
echo "  sbatch slurm/scripts/exp3_throughput.sh"
echo "  sbatch slurm/scripts/exp4_parallel_scaling.sh"
echo "  sbatch slurm/scripts/exp5_ablation.sh"
echo ""
echo "  squeue --me"
echo "======================================================="
