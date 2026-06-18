#!/usr/bin/env bash
# ============================================================
# Download gnomAD-SV v4 structural variant callset
# Used in Experiment 6 to validate biological findings.
#
# Usage:
#   bash experiments/download/download_gnomad_sv.sh
#
# Output: $G2G_DATA_ROOT/gnomad_sv/
# ============================================================
set -euo pipefail

DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
OUT_DIR="$DATA_ROOT/gnomad_sv"
mkdir -p "$OUT_DIR"

echo "[gnomad] Output: $OUT_DIR"

# ─────────────────────────────────────────────────────────────
# gnomAD v4 SV — public GCS bucket (no credentials needed)
# Full SV callset (all chromosomes)
# ─────────────────────────────────────────────────────────────
GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv"

# Main SV VCF (sites only — no genotypes, smaller file ~500 MB)
SV_VCF_URL="${GNOMAD_BASE}/gnomad.v4.1.sv.sites.vcf.gz"
SV_TBI_URL="${SV_VCF_URL}.tbi"

echo "[gnomad] Downloading gnomAD SV v4.1 sites VCF..."
wget -q --show-progress --tries=5 --waitretry=30 \
    -O "$OUT_DIR/gnomad.v4.1.sv.sites.vcf.gz" "$SV_VCF_URL" || {
    # Fallback to v4.0
    SV_VCF_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/genome_sv/gnomad.v4.0.sv.sites.vcf.gz"
    echo "[gnomad] Trying v4.0 fallback..."
    wget -q --show-progress --tries=3 --waitretry=60 \
        -O "$OUT_DIR/gnomad.v4.1.sv.sites.vcf.gz" "$SV_VCF_URL"
}

echo "[gnomad] Downloading tabix index..."
wget -q -O "$OUT_DIR/gnomad.v4.1.sv.sites.vcf.gz.tbi" "$SV_TBI_URL" || \
    tabix -p vcf "$OUT_DIR/gnomad.v4.1.sv.sites.vcf.gz"

# ─────────────────────────────────────────────────────────────
# Extract SV BED file (for bedtools intersection)
# Columns: chrom, start, end, svid, svtype, AF, popmax_af, fst
# ─────────────────────────────────────────────────────────────
echo "[gnomad] Extracting SV BED for intersection..."
bcftools query \
    -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%INFO/AF\t%INFO/POPMAX_AF\n' \
    "$OUT_DIR/gnomad.v4.1.sv.sites.vcf.gz" \
    | awk 'BEGIN{OFS="\t"} {
        # Ensure END >= POS+1
        end = ($3 == ".") ? $2+1 : $3;
        print $1, $2-1, end, $4, $5, $6, $7
    }' \
    | sort -k1,1 -k2,2n \
    > "$OUT_DIR/gnomad.v4.1.sv.bed"

bgzip -@ 4 "$OUT_DIR/gnomad.v4.1.sv.bed"
tabix -p bed "$OUT_DIR/gnomad.v4.1.sv.bed.gz"

echo "[gnomad] SV BED: $OUT_DIR/gnomad.v4.1.sv.bed.gz"

# ─────────────────────────────────────────────────────────────
# Extract high-FST SVs (population-differentiated)
# FST threshold: > 0.15 is conventionally "population-differentiated"
# ─────────────────────────────────────────────────────────────
echo "[gnomad] Extracting high-FST SVs (FST_AFR_EAS > 0.10)..."
bcftools query \
    -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%INFO/AF\t%INFO/POPMAX_AF\t%INFO/FST_AFR_EAS\n' \
    "$OUT_DIR/gnomad.v4.1.sv.sites.vcf.gz" \
    | awk 'BEGIN{OFS="\t"} $8 != "." && $8+0 > 0.10 {
        end = ($3 == ".") ? $2+1 : $3;
        print $1, $2-1, end, $4, $5, $6, $7, $8
    }' \
    | sort -k1,1 -k2,2n \
    > "$OUT_DIR/gnomad.v4.1.sv.high_fst_afr_eas.bed"

echo "[gnomad] High-FST SVs: $(wc -l < "$OUT_DIR/gnomad.v4.1.sv.high_fst_afr_eas.bed") records"

echo "[gnomad] Done. Files in $OUT_DIR:"
ls -lh "$OUT_DIR"
