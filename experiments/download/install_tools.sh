#!/usr/bin/env bash
# ============================================================
# Download / install additional tools not in the base conda env.
# Run after 00_setup_env.sh.
#
#   conda activate g2g_eval
#   bash experiments/download/install_tools.sh
# ============================================================
set -euo pipefail

DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
TOOLS_DIR="$DATA_ROOT/tools"
mkdir -p "$TOOLS_DIR"

echo "[tools] Installing into $TOOLS_DIR"

# ─────────────────────────────────────────────────────────────
# 1. PanGraph (graph-to-graph baseline)
# ─────────────────────────────────────────────────────────────
# NOTE: upstream `neherlab/pangraph` was rewritten in Rust, so the old
# `Pkg.add(url="…/pangraph")` Julia install no longer works — the repo root
# has no Project.toml ("could not find project file"). Bioconda packages a
# ready-to-run `pangraph` with the same `pangraph build …` CLI, which is the
# simplest no-sudo install, so we use that.
echo "[tools] Installing PanGraph via bioconda..."
if conda install -y -c bioconda -c conda-forge pangraph; then
    echo "[tools] PanGraph installed: $(command -v pangraph)"
else
    echo "[tools] WARNING: 'conda install pangraph' failed." >&2
    echo "[tools] Exp3 will still run and record PanGraph rows with exit_code=-1" >&2
    echo "[tools] (skipped). Re-run that baseline later with:" >&2
    echo "[tools]   ONLY_PANGRAPH=1 bash experiments/exp3_throughput/run.sh" >&2
fi

# ─────────────────────────────────────────────────────────────
# 2. minigraph-cactus (for building HPRC pangenomes)
#    Cactus via conda-forge (the only no-sudo method)
# ─────────────────────────────────────────────────────────────
echo "[tools] Installing Cactus (minigraph-cactus pipeline)..."
conda install -y -c bioconda -c conda-forge \
    cactus \
    minigraph \
    gfatools

echo "[tools] Verifying cactus installation..."
cactus --version || echo "[tools] WARNING: cactus not found — may need manual install"

# ─────────────────────────────────────────────────────────────
# 3. gfatools — GFA manipulation (needed for pangenome prep)
# ─────────────────────────────────────────────────────────────
if ! command -v gfatools &>/dev/null; then
    echo "[tools] Building gfatools from source..."
    git clone --depth 1 https://github.com/lh3/gfatools.git "$TOOLS_DIR/gfatools_src"
    make -C "$TOOLS_DIR/gfatools_src"
    cp "$TOOLS_DIR/gfatools_src/gfatools" "$TOOLS_DIR/"
    echo "[tools] gfatools installed at $TOOLS_DIR/gfatools"
fi

# ─────────────────────────────────────────────────────────────
# 4. paftools.js (minimap2 utility — convert PAF to BED etc.)
# ─────────────────────────────────────────────────────────────
if ! command -v paftools.js &>/dev/null; then
    conda install -y -c bioconda -c conda-forge k8 || true
fi

echo ""
echo "[tools] Tool installation complete."
echo "[tools] Tools in PATH: minimap2=$(command -v minimap2 2>/dev/null || echo MISSING)"
echo "[tools]                bedtools=$(command -v bedtools 2>/dev/null || echo MISSING)"
echo "[tools]                pangraph=$(command -v pangraph 2>/dev/null || echo MISSING)"
echo "[tools]                minigraph=$(command -v minigraph 2>/dev/null || echo MISSING)"
