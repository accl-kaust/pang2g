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
# `Pkg.add(url=".../pangraph")` Julia install no longer works (the repo root
# has no Project.toml) and it is not on bioconda. The Rust rewrite ships static
# prebuilt binaries on its GitHub releases, so we fetch the linux x86_64 asset
# directly. The asset URL is discovered from the releases API so a version bump
# doesn't break the install.
echo "[tools] Installing PanGraph (Rust prebuilt binary)..."
install_pangraph() {
    local api="https://api.github.com/repos/neherlab/pangraph/releases/latest"
    local asset_url
    asset_url=$(curl -fsSL "$api" \
        | grep -o '"browser_download_url": *"[^"]*"' \
        | sed 's/.*"browser_download_url": *"//; s/"$//' \
        | grep -iE 'linux' \
        | grep -iE 'x86_64|amd64' \
        | grep -viE '\.sha|\.txt|\.sig|\.asc' \
        | head -1)
    if [[ -z "$asset_url" ]]; then
        echo "[tools] Could not find a linux x86_64 PanGraph asset." >&2
        return 1
    fi
    echo "[tools] PanGraph asset: $asset_url"
    local fname="$TOOLS_DIR/$(basename "$asset_url")"
    curl -fsSL -o "$fname" "$asset_url" || return 1

    case "$fname" in
        *.tar.gz|*.tgz) tar -xzf "$fname" -C "$TOOLS_DIR" ;;
        *.tar.xz)       tar -xJf "$fname" -C "$TOOLS_DIR" ;;
        *.zip)          unzip -o "$fname" -d "$TOOLS_DIR" ;;
        *)              cp "$fname" "$TOOLS_DIR/pangraph" ;;  # raw binary
    esac

    # Expose the executable as $TOOLS_DIR/pangraph
    if [[ ! -x "$TOOLS_DIR/pangraph" ]]; then
        local found
        found=$(find "$TOOLS_DIR" -type f -name 'pangraph*' \
                    ! -name '*.tar.*' ! -name '*.zip' 2>/dev/null | head -1)
        [[ -n "$found" ]] && cp "$found" "$TOOLS_DIR/pangraph"
    fi
    chmod +x "$TOOLS_DIR/pangraph" 2>/dev/null || true
    [[ -x "$TOOLS_DIR/pangraph" ]] || return 1
}

if install_pangraph; then
    export PATH="$TOOLS_DIR:$PATH"
    grep -qsF "export PATH=\"$TOOLS_DIR:\$PATH\"" ~/.bashrc \
        || echo "export PATH=\"$TOOLS_DIR:\$PATH\"" >> ~/.bashrc
    echo "[tools] PanGraph installed: $TOOLS_DIR/pangraph"
    "$TOOLS_DIR/pangraph" --version 2>/dev/null \
        || "$TOOLS_DIR/pangraph" --help 2>/dev/null | head -1 || true
else
    echo "[tools] WARNING: PanGraph install failed." >&2
    echo "[tools] Exp3 still runs and records PanGraph rows with exit_code=-1" >&2
    echo "[tools] (skipped). Backfill later once installed with:" >&2
    echo "[tools]   ONLY_PANGRAPH=1 bash experiments/exp3_throughput/run.sh" >&2
fi

# ─────────────────────────────────────────────────────────────
# 2. minigraph-cactus (for building HPRC pangenomes)
#    Cactus via conda-forge (the only no-sudo method)
# ─────────────────────────────────────────────────────────────
echo "[tools] Installing Cactus (minigraph-cactus pipeline)..."
# mash is a runtime dependency of cactus-pangenome's minigraph step
# (cactus shells out to `mash sketch`); the bioconda cactus package does NOT
# pull it in, so install it explicitly or the build dies with exit 127.
conda install -y -c bioconda -c conda-forge \
    cactus \
    minigraph \
    gfatools \
    mash

echo "[tools] Verifying cactus installation..."
cactus --version || echo "[tools] WARNING: cactus not found — may need manual install"
for t in mash minigraph gfatools; do
    command -v "$t" &>/dev/null && echo "[tools] $t: $(command -v $t)" || \
        echo "[tools] WARNING: $t not found on PATH"
done

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
