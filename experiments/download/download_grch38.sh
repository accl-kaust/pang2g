#!/usr/bin/env bash
# ============================================================
# Download the GRCh38 reference (no-alt analysis set) used as the
# --reference for cactus-pangenome in Experiment 6 (stage 02_build).
#
# Run on the LOGIN NODE (compute nodes have no internet).
#
# Usage:
#   bash experiments/download/download_grch38.sh
#
# Output: $G2G_DATA_ROOT/reference/GRCh38.fa.gz
# ============================================================
set -euo pipefail

DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
OUT_DIR="$DATA_ROOT/reference"
OUT="$OUT_DIR/GRCh38.fa.gz"
mkdir -p "$OUT_DIR"

# The HPRC pangenome pipeline uses the GRCh38 "no-alt analysis set".
URLS=(
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
)

echo "[grch38] Output: $OUT"
if [[ -f "$OUT" ]] && [[ $(stat -c%s "$OUT") -gt 700000000 ]]; then
    echo "[grch38] SKIP download (reference already present, $(stat -c%s "$OUT") bytes)."
else
    ok=0
    for u in "${URLS[@]}"; do
        echo "[grch38] Trying: $u"
        if wget -q --show-progress --tries=3 --waitretry=20 -O "$OUT.tmp" "$u"; then
            # Ensure the result is bgzip/gzip-compressed; if a plain .fa was
            # fetched, compress it so downstream tools get a .fa.gz.
            if file "$OUT.tmp" | grep -qi 'gzip'; then
                mv "$OUT.tmp" "$OUT"
            else
                echo "[grch38] Compressing plain FASTA -> bgzip..."
                bgzip -c "$OUT.tmp" > "$OUT" && rm -f "$OUT.tmp"
            fi
            ok=1; break
        fi
    done
    if [[ "$ok" -eq 0 ]]; then
        echo "[grch38] ERROR: could not download GRCh38 from any mirror." >&2
        rm -f "$OUT.tmp"
        exit 1
    fi
fi

echo "[grch38] Indexing (samtools faidx)..."
samtools faidx "$OUT" 2>/dev/null || \
    echo "[grch38] WARNING: faidx failed (a bgzipped FASTA is required for faidx)."

echo "[grch38] Done:"
ls -lh "$OUT"*
