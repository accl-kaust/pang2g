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
# 1. PanGraph (Julia package — bacterial graph-to-graph baseline)
# ─────────────────────────────────────────────────────────────
echo "[tools] Installing PanGraph via Julia..."

# Install Julia if not already on PATH
if ! command -v julia &>/dev/null; then
    echo "[tools] Julia not found — downloading Julia 1.10..."
    JULIA_TGZ="$TOOLS_DIR/julia.tar.gz"
    JULIA_URL="https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.4-linux-x86_64.tar.gz"
    wget -q --show-progress -O "$JULIA_TGZ" "$JULIA_URL"
    tar -xzf "$JULIA_TGZ" -C "$TOOLS_DIR"
    JULIA_BIN=$(ls -d "$TOOLS_DIR"/julia-*/bin | head -1)
    export PATH="$JULIA_BIN:$PATH"
    echo "export PATH=\"$JULIA_BIN:\$PATH\"" >> ~/.bashrc
    echo "[tools] Julia installed at $JULIA_BIN"
fi

# Install PanGraph
julia -e '
using Pkg
Pkg.add(url="https://github.com/neherlab/pangraph")
Pkg.build("PanGraph")
'
# Create a wrapper script
cat > "$TOOLS_DIR/pangraph" <<'WRAP'
#!/usr/bin/env bash
julia -e "using PanGraph; PanGraph.main()" -- "$@"
WRAP
chmod +x "$TOOLS_DIR/pangraph"
export PATH="$TOOLS_DIR:$PATH"
echo "export PATH=\"$TOOLS_DIR:\$PATH\"" >> ~/.bashrc
echo "[tools] PanGraph wrapper at $TOOLS_DIR/pangraph"

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
echo "[tools]                julia=$(command -v julia 2>/dev/null || echo MISSING)"
echo "[tools]                minigraph=$(command -v minigraph 2>/dev/null || echo MISSING)"
