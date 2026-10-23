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

# Fail loudly if the toolchain isn't on PATH (usual cause: forgot to
# `conda activate g2g_eval`) rather than silently producing empty output.
for tool in bcftools bgzip tabix wget; do
    if ! command -v "$tool" &>/dev/null; then
        echo "[gnomad] ERROR: '$tool' not found on PATH." >&2
        echo "[gnomad]        Did you 'conda activate g2g_eval' first?" >&2
        exit 1
    fi
done

echo "[gnomad] Output: $OUT_DIR"

# -------------------------------------------------------------
# gnomAD v4 SV - public GCS bucket (no credentials needed)
# Full SV callset (all chromosomes)
# -------------------------------------------------------------
GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv"

# Main SV VCF (sites only - no genotypes, smaller file ~500 MB)
SV_VCF_URL="${GNOMAD_BASE}/gnomad.v4.1.sv.sites.vcf.gz"
SV_TBI_URL="${SV_VCF_URL}.tbi"

if [[ -f "$OUT_DIR/gnomad.v4.1.sv.sites.vcf.gz" ]] && \
   [[ $(stat -c%s "$OUT_DIR/gnomad.v4.1.sv.sites.vcf.gz") -gt 100000000 ]]; then
    echo "[gnomad] SKIP download (VCF already present)."
else
    echo "[gnomad] Downloading gnomAD SV v4.1 sites VCF..."
    wget -q --show-progress --tries=5 --waitretry=30 \
        -O "$OUT_DIR/gnomad.v4.1.sv.sites.vcf.gz" "$SV_VCF_URL" || {
        # Fallback to v4.0
        SV_VCF_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/genome_sv/gnomad.v4.0.sv.sites.vcf.gz"
        echo "[gnomad] Trying v4.0 fallback..."
        wget -q --show-progress --tries=3 --waitretry=60 \
            -O "$OUT_DIR/gnomad.v4.1.sv.sites.vcf.gz" "$SV_VCF_URL"
    }
fi

echo "[gnomad] Downloading tabix index..."
wget -q -O "$OUT_DIR/gnomad.v4.1.sv.sites.vcf.gz.tbi" "$SV_TBI_URL" || \
    tabix -p vcf "$OUT_DIR/gnomad.v4.1.sv.sites.vcf.gz"

# -------------------------------------------------------------
# Discover which INFO tags this release actually carries.
# gnomAD v4.1 SV does NOT have INFO/POPMAX_AF (popmax was renamed to grpmax
# in v4) and ships no per-pair FST field at all - so we detect the real tag
# names from the header and compute AFR/EAS FST ourselves.
# -------------------------------------------------------------
VCF="$OUT_DIR/gnomad.v4.1.sv.sites.vcf.gz"
HDR="$(bcftools view -h "$VCF")"
have_tag() { grep -q "##INFO=<ID=$1," <<< "$HDR"; }

# popmax/grpmax AF column (optional)
PMAX_TAG=""
for t in grpmax_AF AF_grpmax POPMAX_AF AF_popmax grpmax popmax; do
    have_tag "$t" && { PMAX_TAG="$t"; break; }
done
# per-population AF tags for the AFR/EAS FST (try common casings)
AFR_TAG=""; EAS_TAG=""
for t in AF_afr AF_AFR afr_AF AFR_AF; do have_tag "$t" && { AFR_TAG="$t"; break; }; done
for t in AF_eas AF_EAS eas_AF EAS_AF; do have_tag "$t" && { EAS_TAG="$t"; break; }; done
echo "[gnomad] tags: popmax='${PMAX_TAG:-none}' afr='${AFR_TAG:-none}' eas='${EAS_TAG:-none}'"

# -------------------------------------------------------------
# Extract SV BED file (for bedtools intersection)
# Columns: chrom, start, end, svid, svtype, AF, popmax_af(. if absent)
# -------------------------------------------------------------
echo "[gnomad] Extracting SV BED for intersection..."
PMAX_FMT='.'; [[ -n "$PMAX_TAG" ]] && PMAX_FMT="%INFO/$PMAX_TAG"
bcftools query \
    -f "%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%INFO/AF\t${PMAX_FMT}\n" \
    "$VCF" \
    | awk 'BEGIN{OFS="\t"} {
        end = ($3 == ".") ? $2+1 : $3;
        print $1, $2-1, end, $4, $5, $6, $7
    }' \
    | sort -k1,1 -k2,2n \
    > "$OUT_DIR/gnomad.v4.1.sv.bed"

bgzip -f -@ 4 "$OUT_DIR/gnomad.v4.1.sv.bed"
tabix -f -p bed "$OUT_DIR/gnomad.v4.1.sv.bed.gz"

echo "[gnomad] SV BED: $OUT_DIR/gnomad.v4.1.sv.bed.gz"

# -------------------------------------------------------------
# Extract population-differentiated SVs (high AFR/EAS FST).
# gnomAD ships no FST field, so compute Wright's FST from the per-population
# allele frequencies: pbar=(p_afr+p_eas)/2; Ht=2*pbar*(1-pbar);
# Hs=(2*p_afr*(1-p_afr)+2*p_eas*(1-p_eas))/2; FST=(Ht-Hs)/Ht.
# Threshold > 0.10 ≈ population-differentiated.
# -------------------------------------------------------------
HIFST_BED="$OUT_DIR/gnomad.v4.1.sv.high_fst_afr_eas.bed"
if [[ -n "$AFR_TAG" && -n "$EAS_TAG" ]]; then
    echo "[gnomad] Computing AFR/EAS FST from $AFR_TAG, $EAS_TAG (threshold > 0.10)..."
    bcftools query \
        -f "%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%INFO/AF\t%INFO/$AFR_TAG\t%INFO/$EAS_TAG\n" \
        "$VCF" \
        | awk 'BEGIN{OFS="\t"}
            $7 != "." && $8 != "." {
                a=$7+0; e=$8+0;
                pbar=(a+e)/2; Ht=2*pbar*(1-pbar);
                if (Ht <= 0) next;
                Hs=(2*a*(1-a)+2*e*(1-e))/2;
                fst=(Ht-Hs)/Ht;
                if (fst > 0.10) {
                    end = ($3 == ".") ? $2+1 : $3;
                    printf "%s\t%d\t%s\t%s\t%s\t%s\t%.4f\n", $1, $2-1, end, $4, $5, $6, fst
                }
            }' \
        | sort -k1,1 -k2,2n \
        > "$HIFST_BED"
    echo "[gnomad] High-FST SVs: $(wc -l < "$HIFST_BED") records"
else
    echo "[gnomad] WARNING: per-population AF tags not found - skipping high-FST BED." >&2
    : > "$HIFST_BED"
fi

echo "[gnomad] Done. Files in $OUT_DIR:"
ls -lh "$OUT_DIR"
