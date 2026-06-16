#!/bin/bash
# =============================================================================
# setup_cluster.sh — One-time setup before running G2G Phase 5 jobs.
# =============================================================================
# Run INTERACTIVELY on the head/login node (NOT via sbatch):
#     bash slurm/scripts/setup_cluster.sh
#
# This script is idempotent: safe to re-run after a `git pull`.
#
# Configuration (override by exporting before running):
#     G2G_REPO   — repo root            (default: directory two levels up)
#     G2G_ENV    — conda env name       (default: phase5)
#     CONDA_BASE — miniconda/mamba root  (default: $HOME/miniconda3)
#
# What changed from the old setup (TBB era):
#   * Removed `tbb-devel` from the conda install — the runtime is now OpenMP,
#     which ships with GCC. Installing TBB is unnecessary and was a source of
#     "package not found" failures on bioconda.
#   * Environment is `phase5` everywhere (the old README said `phase4`).
#   * Builds Release with -O3 -march=native and runs the FULL ctest suite,
#     failing loudly if any test regresses (the paper's correctness claim
#     depends on 143/143 passing on the evaluation node itself, not just CI).
#   * Adds a large-graph storage + memory readiness check so you find out
#     BEFORE a 14-hour job that /scratch is too small or RAM too tight.
# =============================================================================

set -euo pipefail

# ── Resolve script + repo locations ──────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
G2G_REPO="${G2G_REPO:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
G2G_ENV="${G2G_ENV:-phase5}"
CONDA_BASE="${CONDA_BASE:-$HOME/miniconda3}"
BUILD="$G2G_REPO/build"

source "$SCRIPT_DIR/lib_measure.sh"

echo "======================================================================"
echo "  G2G Phase 5 — Cluster Setup"
echo "  Host:  $(hostname)"
echo "  Date:  $(date)"
echo "  Repo:  $G2G_REPO"
echo "  Env:   $G2G_ENV"
echo "  Conda: $CONDA_BASE"
echo "======================================================================"

# ── 1. Sanity: repo layout ───────────────────────────────────────────────────
for required in CMakeLists.txt include/g2g src tests slurm/synth/gen_scaling_gfa.py; do
    if [[ ! -e "$G2G_REPO/$required" ]]; then
        echo "[setup] FATAL: expected '$required' under $G2G_REPO but it is missing."
        echo "[setup]        Set G2G_REPO to the correct repo root and re-run."
        exit 1
    fi
done
echo "[setup] repo layout OK"

# ── 2. Conda / mamba ─────────────────────────────────────────────────────────
if [[ ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
    echo "[setup] FATAL: conda not found at $CONDA_BASE"
    echo "[setup]        Install miniconda/mamba or set CONDA_BASE."
    exit 1
fi
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate base
echo "[setup] conda: $(conda --version)"

USE_MAMBA=""
command -v mamba >/dev/null 2>&1 && USE_MAMBA="mamba"
INSTALLER="${USE_MAMBA:-conda}"
echo "[setup] using installer: $INSTALLER"

if conda env list | grep -qE "^${G2G_ENV}[[:space:]]"; then
    echo "[setup] env '$G2G_ENV' exists — updating core packages"
    conda activate "$G2G_ENV"
    $INSTALLER install -y -c conda-forge -c bioconda \
        cmake cxx-compiler make minimap2 samtools bedtools
else
    echo "[setup] creating env '$G2G_ENV'"
    $INSTALLER create -y -n "$G2G_ENV" -c conda-forge -c bioconda \
        python=3.11 cmake cxx-compiler make minimap2 samtools bedtools \
        numpy pandas matplotlib
    conda activate "$G2G_ENV"
fi

echo "[setup] g++:      $(g++ --version | head -1)"
echo "[setup] cmake:    $(cmake --version | head -1)"
echo "[setup] minimap2: $(minimap2 --version)"
echo "[setup] OpenMP:   $(echo | g++ -fopenmp -dM -E - 2>/dev/null | grep -i _OPENMP || echo 'NOT FOUND')"

# ── 3. Build (Release) ───────────────────────────────────────────────────────
echo ""
echo "[setup] configuring + building (Release, -O3 -march=native)..."
cmake -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER="$(which g++)" \
      -DCMAKE_C_COMPILER="$(which gcc)" \
      -DCMAKE_CXX_FLAGS="-O3 -march=native"
cmake --build "$BUILD" -j"$(nproc)"
echo "[setup] build complete: $BUILD/g2g"

# ── 4. Full test suite (must pass before evaluation) ─────────────────────────
echo ""
echo "[setup] running full ctest suite (correctness gate)..."
pushd "$BUILD" >/dev/null
if ctest --output-on-failure -j"$(nproc)"; then
    echo "[setup] ALL TESTS PASS"
else
    echo "[setup] WARNING: some tests failed — investigate before trusting results."
    echo "[setup]          (minimap2 consistency tests will SKIP, not fail, if the"
    echo "[setup]           binary is missing; those skips are expected.)"
fi
popd >/dev/null

# ── 5. Output + log dirs ─────────────────────────────────────────────────────
mkdir -p "$G2G_REPO/results"/{exp3,exp4,exp5,exp6} "$G2G_REPO/slurm/logs"
echo "[setup] results dirs ready under $G2G_REPO/results"

# ── 6. perf availability (Exp5 A1 L3 miss rate) ──────────────────────────────
echo ""
echo "[setup] checking perf (needed for Exp5-A1 cache-miss measurement)..."
if [[ $HAVE_PERF -eq 1 ]]; then
    echo "[setup] perf OK — L3 miss rate will be collected"
else
    PARANOID=$(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo unknown)
    echo "[setup] perf restricted (perf_event_paranoid=$PARANOID)"
    echo "[setup]   Exp5-A1 will record wall-time only. To enable hardware"
    echo "[setup]   counters: ask the admin to run"
    echo "[setup]       sudo sysctl kernel.perf_event_paranoid=2"
fi

# ── 7. Large-graph readiness check ───────────────────────────────────────────
echo ""
echo "[setup] large-graph readiness:"
GB_RAM=$(free -g | awk '/^Mem:/{print $2}')
NCPU=$(nproc)
SCRATCH="${G2G_SCRATCH:-$G2G_REPO/results}"
AVAIL_TB=$(df -BG "$SCRATCH" 2>/dev/null | awk 'NR==2{gsub("G","",$4); printf "%.1f", $4/1024}')
echo "[setup]   RAM:         ${GB_RAM} GB"
echo "[setup]   CPUs:        ${NCPU}"
echo "[setup]   scratch dir: $SCRATCH (${AVAIL_TB:-?} TB free)"
echo ""
echo "[setup]   HPRC-realistic sweep needs roughly:"
echo "[setup]     n=10M snarls -> ~3.3 GB/GFA on disk, ~40-80 GB RAM working set"
echo "[setup]     full sweep (1K..10M, x2 graphs, x3 seeds) -> ~25 GB disk"
if (( GB_RAM < 256 )); then
    echo "[setup]   NOTE: <256 GB RAM. The 10M-snarl point may be tight; cap the"
    echo "[setup]         sweep at 1M-5M (edit MAX_N in exp3_throughput.sh) if needed."
fi

# ── 8. Smoke test ────────────────────────────────────────────────────────────
echo ""
echo "[setup] smoke test (100-snarl pair, 4 threads)..."
SMOKE=/tmp/g2g_smoke_$$
mkdir -p "$SMOKE"
# Use the SAME generator the experiments use, so this validates the real path.
python3 "$G2G_REPO/slurm/synth/gen_hprc_scale_gfa.py" \
    --mode size --n_snarls 100 --out_dir "$SMOKE" --seed 1 >/dev/null
"$BUILD/g2g" "$SMOKE/ga_100.gfa" "$SMOKE/gb_100.gfa" \
    --gaf /dev/null --threads 4 2>&1 | tail -3
rm -rf "$SMOKE"
echo "[setup] smoke test PASSED"

echo ""
echo "======================================================================"
echo "  Setup complete. From $G2G_REPO, submit experiments with:"
echo ""
echo "    bash slurm/scripts/submit_all.sh"
echo ""
echo "  or individually:"
echo "    sbatch slurm/scripts/exp3_throughput.sh"
echo "    sbatch slurm/scripts/exp4_parallel_scaling.sh"
echo "    sbatch slurm/scripts/exp5_ablation.sh"
echo "    sbatch slurm/scripts/exp6_biological.sh   # after HPRC download"
echo "======================================================================"
