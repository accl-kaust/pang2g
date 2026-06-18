#!/usr/bin/env bash
# ============================================================
# Phase 5 — Environment Setup
# Creates the g2g_eval conda environment with all tools needed
# for all six experiments.
#
# Run once before any experiment:
#   bash experiments/00_setup_env.sh
#
# No sudo required — uses miniconda.
# ============================================================
set -euo pipefail

ENV_NAME="g2g_eval"
G2G_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "[setup] G2G root: $G2G_ROOT"
echo "[setup] Creating conda env: $ENV_NAME"

# ── Create base env ───────────────────────────────────────────
conda create -y -n "$ENV_NAME" python=3.11

# ── Core build tools ──────────────────────────────────────────
conda run -n "$ENV_NAME" conda install -y -c conda-forge \
    cmake=3.27 \
    make \
    ninja \
    cxx-compiler \
    c-compiler \
    openmp \
    gdb \
    valgrind

# ── Bioinformatics tools ──────────────────────────────────────
conda run -n "$ENV_NAME" conda install -y -c bioconda -c conda-forge \
    minimap2 \
    bedtools \
    samtools \
    bcftools \
    seqkit \
    paftools

# ── Python analysis packages ──────────────────────────────────
conda run -n "$ENV_NAME" pip install \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    scipy \
    networkx \
    memory_profiler \
    psutil \
    tqdm \
    pysam

# ── perf stat wrapper (user-space; no sudo needed for basic events) ───
# perf is in the kernel tools, should be available on HPC nodes
# We alias it to use --per-thread for user-space profiling
cat >> ~/.bashrc <<'EOF'
export PATH="$CONDA_PREFIX/bin:$PATH"
EOF

echo ""
echo "[setup] Environment '$ENV_NAME' ready."
echo "[setup] Activate with:  conda activate $ENV_NAME"
echo "[setup] Then build with: bash experiments/01_build_all.sh"
