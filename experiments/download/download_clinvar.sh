#!/usr/bin/env bash
# ============================================================
# Download ClinVar VCF (pathogenic/likely pathogenic SVs)
# Used in Experiment 6 to check if shared G2G alleles contain
# known pathogenic variants present across populations.
#
# Usage:
#   bash experiments/download/download_clinvar.sh
#
# Output: $G2G_DATA_ROOT/clinvar/
# ============================================================
set -euo pipefail

DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
OUT_DIR="$DATA_ROOT/clinvar"
mkdir -p "$OUT_DIR"

NCBI_FTP="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"

echo "[clinvar] Output: $OUT_DIR"

# ─────────────────────────────────────────────────────────────
# 1. Download ClinVar VCF (GRCh38)
# ─────────────────────────────────────────────────────────────
echo "[clinvar] Downloading ClinVar VCF (GRCh38)..."
wget -q --show-progress --tries=5 --waitretry=30 \
    -O "$OUT_DIR/clinvar.vcf.gz" \
    "${NCBI_FTP}/clinvar.vcf.gz"

wget -q -O "$OUT_DIR/clinvar.vcf.gz.tbi" \
    "${NCBI_FTP}/clinvar.vcf.gz.tbi" || \
    tabix -p vcf "$OUT_DIR/clinvar.vcf.gz"

# ─────────────────────────────────────────────────────────────
# 2. Extract pathogenic/likely pathogenic SVs
#    CLNSIG = Pathogenic or Likely_pathogenic
#    CLNVC  = structural variant types (Deletion, Duplication, etc.)
# ─────────────────────────────────────────────────────────────
echo "[clinvar] Extracting pathogenic SVs..."
bcftools query \
    -i 'CLNSIG~"Pathogenic" && (CLNVC~"deletion" || CLNVC~"duplication" || CLNVC~"inversion" || CLNVC~"insertion")' \
    -f '%CHROM\t%POS\t%END\t%ID\t%CLNVC\t%CLNSIG\t%CLNDN\n' \
    "$OUT_DIR/clinvar.vcf.gz" \
    | awk 'BEGIN{OFS="\t"} {
        end = ($3 == ".") ? $2+1 : $3;
        print $1, $2-1, end, $4, $5, $6, $7
    }' \
    | sort -k1,1 -k2,2n \
    > "$OUT_DIR/clinvar_pathogenic_svs.bed"

# ─────────────────────────────────────────────────────────────
# 3. All ClinVar SVs (any significance) — for sensitivity analysis
# ─────────────────────────────────────────────────────────────
echo "[clinvar] Extracting all ClinVar SVs..."
bcftools query \
    -i '(CLNVC~"deletion" || CLNVC~"duplication" || CLNVC~"inversion" || CLNVC~"insertion")' \
    -f '%CHROM\t%POS\t%END\t%ID\t%CLNVC\t%CLNSIG\t%CLNDN\n' \
    "$OUT_DIR/clinvar.vcf.gz" \
    | awk 'BEGIN{OFS="\t"} {
        end = ($3 == ".") ? $2+1 : $3;
        print $1, $2-1, end, $4, $5, $6, $7
    }' \
    | sort -k1,1 -k2,2n \
    > "$OUT_DIR/clinvar_all_svs.bed"

# ─────────────────────────────────────────────────────────────
# 4. Create GRCh38 chromosome size file (needed for bedtools)
# ─────────────────────────────────────────────────────────────
echo "[clinvar] Creating chromosome sizes file..."
wget -q -O "$OUT_DIR/hg38.chrom.sizes" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes" || \
    echo "[clinvar] WARNING: Could not download chrom sizes (needed for bedtools slop)"

echo ""
echo "[clinvar] Summary:"
echo "  Pathogenic SVs: $(wc -l < "$OUT_DIR/clinvar_pathogenic_svs.bed")"
echo "  All SVs:        $(wc -l < "$OUT_DIR/clinvar_all_svs.bed")"
echo "[clinvar] Files in $OUT_DIR:"
ls -lh "$OUT_DIR"
