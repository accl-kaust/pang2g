#!/usr/bin/env bash
# Download the third-party tools used by the REAL-data demo (04) into demo/tools/.
# The synthetic demo (03) needs none of these - only g2g + python3.
#
# Tools:
#   vg        - build pangenome graphs from a reference + VCF (vgteam)
#   minimap2  - align haplotypes to the reference (lh3)
#   paftools.js + k8 - call variants from the alignment (ships with minimap2)
#   pysam     - bgzip/tabix/faidx (pip)
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

VG_VER=v1.75.0
MM2_TAG=v2.28
MM2_FILE=minimap2-2.28_x64-linux   # release tag carries a 'v', the asset does not

log "Installing tools into $TOOLS"

if [ ! -x "$TOOLS/vg" ]; then
    note "downloading vg $VG_VER (~52 MB)…"
    curl -fsSL --retry 3 -o "$TOOLS/vg" \
        "https://github.com/vgteam/vg/releases/download/$VG_VER/vg"
    chmod +x "$TOOLS/vg"
else note "vg already present"; fi

if [ ! -x "$TOOLS/minimap2" ]; then
    note "downloading minimap2 $MM2_TAG (+ paftools.js, k8)…"
    tmp="$WORK/mm2.tar.bz2"
    curl -fsSL --retry 3 -o "$tmp" \
        "https://github.com/lh3/minimap2/releases/download/$MM2_TAG/$MM2_FILE.tar.bz2"
    tar xjf "$tmp" -C "$WORK"
    cp "$WORK/$MM2_FILE/"{minimap2,k8,paftools.js} "$TOOLS/"
else note "minimap2 already present"; fi

note "ensuring pysam (bgzip/tabix/faidx) is importable…"
python3 -c "import pysam" 2>/dev/null || pip install --quiet pysam

log "Tool versions"
"$TOOLS/vg" version | head -1
"$TOOLS/minimap2" --version | sed 's/^/minimap2 /'
python3 -c "import pysam; print('pysam', pysam.__version__)"
echo "OK - tools ready in $TOOLS"
