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

# Fail loudly if the toolchain isn't on PATH — otherwise bcftools errors get
# swallowed by `|| true` and we silently write empty BEDs. The usual cause is
# forgetting to `conda activate g2g_eval` before running.
for tool in bcftools tabix wget; do
    if ! command -v "$tool" &>/dev/null; then
        echo "[clinvar] ERROR: '$tool' not found on PATH." >&2
        echo "[clinvar]        Did you 'conda activate g2g_eval' first?" >&2
        exit 1
    fi
done

echo "[clinvar] Output: $OUT_DIR"

# ─────────────────────────────────────────────────────────────
# 1. Download ClinVar VCF (GRCh38)
# ─────────────────────────────────────────────────────────────
echo "[clinvar] Downloading ClinVar VCF (GRCh38)..."
if [[ -f "$OUT_DIR/clinvar.vcf.gz" ]] && \
   [[ $(stat -c%s "$OUT_DIR/clinvar.vcf.gz") -gt 50000000 ]] && \
   bcftools view -h "$OUT_DIR/clinvar.vcf.gz" &>/dev/null; then
    echo "[clinvar] SKIP download (valid VCF already present, $(stat -c%s "$OUT_DIR/clinvar.vcf.gz") bytes)."
else
    wget -q --show-progress --tries=5 --waitretry=30 \
        -O "$OUT_DIR/clinvar.vcf.gz" \
        "${NCBI_FTP}/clinvar.vcf.gz"
fi

wget -q -O "$OUT_DIR/clinvar.vcf.gz.tbi" \
    "${NCBI_FTP}/clinvar.vcf.gz.tbi" || \
    tabix -p vcf "$OUT_DIR/clinvar.vcf.gz"

# ─────────────────────────────────────────────────────────────
# 2. Extract pathogenic/likely pathogenic SVs
#    CLNSIG = Pathogenic or Likely_pathogenic
#    CLNVC  = structural variant types (Deletion, Duplication, etc.)
# ─────────────────────────────────────────────────────────────
# Show top CLNVC values to confirm filter will match
echo "[clinvar] CLNVC value counts (top 10):"
bcftools query -f '%CLNVC\n' "$OUT_DIR/clinvar.vcf.gz" 2>/dev/null \
    | sort | uniq -c | sort -rn | head -10 || true

echo "[clinvar] Extracting pathogenic SVs..."
# CLNVC uses Title_Case (e.g. "Deletion", "Duplication"); bcftools ~ is case-sensitive.
# Use case-insensitive regex via (?i) prefix (POSIX ERE supported by bcftools >= 1.10).
# Do NOT filter on FILTER column — ClinVar variants routinely carry filters like
# "no_assertion_criteria_provided" but are still clinically relevant.
SV_FILTER='CLNVC~"(?i)(deletion|duplication|inversion|insertion|copy_number)"'
PATHO_FILTER='CLNSIG~"(?i)(pathogenic)"'
bcftools query \
    -i "${PATHO_FILTER} && ${SV_FILTER}" \
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
    -i "${SV_FILTER}" \
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
# Try multiple sources; fall back to generating from the fai of the GRCh38 reference.
SIZES_OK=0
for url in \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes" \
    "https://genome.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes" \
    "https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dict"
do
    if wget -q --tries=3 --waitretry=10 -O "$OUT_DIR/hg38.chrom.sizes.tmp" "$url" 2>/dev/null && \
       [[ $(stat -c%s "$OUT_DIR/hg38.chrom.sizes.tmp" 2>/dev/null || echo 0) -gt 1000 ]]; then
        # The .dict file needs to be converted to chrom.sizes format
        if [[ "$url" == *.dict ]]; then
            grep '^@SQ' "$OUT_DIR/hg38.chrom.sizes.tmp" \
                | sed 's/.*SN:\([^\t]*\).*LN:\([^\t]*\).*/\1\t\2/' \
                > "$OUT_DIR/hg38.chrom.sizes"
        else
            mv "$OUT_DIR/hg38.chrom.sizes.tmp" "$OUT_DIR/hg38.chrom.sizes"
        fi
        rm -f "$OUT_DIR/hg38.chrom.sizes.tmp"
        SIZES_OK=1
        echo "[clinvar] chrom.sizes downloaded ($(wc -l < "$OUT_DIR/hg38.chrom.sizes") contigs)"
        break
    fi
done
rm -f "$OUT_DIR/hg38.chrom.sizes.tmp"
if [[ "$SIZES_OK" -eq 0 ]]; then
    # Last resort: extract from the ClinVar VCF header contig lines
    bcftools view -h "$OUT_DIR/clinvar.vcf.gz" \
        | grep '^##contig' \
        | sed 's/.*ID=\([^,>]*\).*length=\([0-9]*\).*/\1\t\2/' \
        | grep -v '^$' \
        > "$OUT_DIR/hg38.chrom.sizes" || true
    if [[ $(stat -c%s "$OUT_DIR/hg38.chrom.sizes" 2>/dev/null || echo 0) -gt 0 ]]; then
        echo "[clinvar] chrom.sizes extracted from VCF header ($(wc -l < "$OUT_DIR/hg38.chrom.sizes") contigs)"
        SIZES_OK=1
    else
        echo "[clinvar] WARNING: could not create hg38.chrom.sizes — bedtools slop will fail" >&2
    fi
fi

echo ""
echo "[clinvar] Summary:"
echo "  Pathogenic SVs: $(wc -l < "$OUT_DIR/clinvar_pathogenic_svs.bed")"
echo "  All SVs:        $(wc -l < "$OUT_DIR/clinvar_all_svs.bed")"
echo "[clinvar] Files in $OUT_DIR:"
ls -lh "$OUT_DIR"
