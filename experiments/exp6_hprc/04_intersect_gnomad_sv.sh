#!/usr/bin/env bash
# ============================================================
# 04_intersect_gnomad_sv.sh - Intersect G2G delta alleles with gnomAD-SV
#
# Tests whether G2G-detected population-specific alleles (delta
# graph) are enriched for population-differentiated SVs in gnomAD.
#
# Hypothesis:
#   AFR-only delta alleles should overlap gnomAD SVs with high
#   AFR–EAS FST (population differentiation).
#
# Produces:
#   - Overlap BED (which delta alleles hit gnomAD SVs)
#   - Enrichment table (observed vs. expected overlap rate)
#   - FST distribution of overlapping vs. non-overlapping delta alleles
#
# Usage:
#   bash experiments/exp6_hprc/04_intersect_gnomad_sv.sh
# ============================================================
set -euo pipefail

DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
RESULTS="$DATA_ROOT/results/exp6"
GNOMAD_DIR="$DATA_ROOT/gnomad_sv"
TIMESTAMP=$(date +%Y%m%dT%H%M%S)

mkdir -p "$RESULTS/intersections"

echo "[gnomad_int] Results: $RESULTS/intersections"

# -------------------------------------------------------------
# Find the most recent G2G delta BED
# -------------------------------------------------------------
DELTA_BED=$(ls -t "$RESULTS"/afr_vs_eas_delta_*_sorted.bed.gz 2>/dev/null | head -1)
if [[ -z "$DELTA_BED" ]]; then
    echo "[gnomad_int] ERROR: No delta BED found in $RESULTS"
    echo "[gnomad_int] Run 03_run_g2g_on_hprc.sh first."
    exit 1
fi
echo "[gnomad_int] Delta BED: $DELTA_BED"

# gnomAD SV files
GNOMAD_BED="$GNOMAD_DIR/gnomad.v4.1.sv.bed.gz"
GNOMAD_HIFST="$GNOMAD_DIR/gnomad.v4.1.sv.high_fst_afr_eas.bed"

if [[ ! -f "$GNOMAD_BED" ]]; then
    echo "[gnomad_int] ERROR: gnomAD BED not found: $GNOMAD_BED"
    echo "[gnomad_int] Run download/download_gnomad_sv.sh first."
    exit 1
fi

# -------------------------------------------------------------
# 1. Basic overlap: delta alleles that hit any gnomAD SV
# -------------------------------------------------------------
OUT_OVERLAP="$RESULTS/intersections/delta_vs_gnomad_${TIMESTAMP}.bed"
echo "[gnomad_int] Computing overlap with all gnomAD SVs..."

bedtools intersect \
    -a "$DELTA_BED" \
    -b "$GNOMAD_BED" \
    -wa -wb \
    > "$OUT_OVERLAP"

N_DELTA=$(zcat "$DELTA_BED" | wc -l)
N_OVERLAP=$(wc -l < "$OUT_OVERLAP")
echo "[gnomad_int] Delta alleles:      $N_DELTA"
echo "[gnomad_int] Overlapping gnomAD: $N_OVERLAP ($(python3 -c "print(f'{$N_OVERLAP/$N_DELTA*100:.1f}%')" 2>/dev/null || echo "?"))"

# -------------------------------------------------------------
# 2. Enrichment for high-FST (AFR-EAS differentiated) SVs
# -------------------------------------------------------------
OUT_HIFST_OVERLAP="$RESULTS/intersections/delta_vs_gnomad_hifst_${TIMESTAMP}.bed"
echo "[gnomad_int] Computing overlap with high-FST gnomAD SVs..."

bedtools intersect \
    -a "$DELTA_BED" \
    -b "$GNOMAD_HIFST" \
    -wa -wb \
    > "$OUT_HIFST_OVERLAP"

N_HIFST_OVERLAP=$(wc -l < "$OUT_HIFST_OVERLAP")
echo "[gnomad_int] High-FST overlap:   $N_HIFST_OVERLAP"

# -------------------------------------------------------------
# 3. Background overlap rate (random BED regions, same size)
# -------------------------------------------------------------
echo "[gnomad_int] Computing background overlap rate..."
CHR_SIZES="$DATA_ROOT/clinvar/hg38.chrom.sizes"

if [[ -f "$CHR_SIZES" ]]; then
    SHUFFLED="$RESULTS/intersections/delta_shuffled_${TIMESTAMP}.bed"
    bedtools shuffle \
        -i "$DELTA_BED" \
        -g "$CHR_SIZES" \
        -noOverlapping \
        > "$SHUFFLED" 2>/dev/null

    SHUFFLED_OVERLAP=$(bedtools intersect -a "$SHUFFLED" -b "$GNOMAD_BED" -u | wc -l || echo 0)
    N_SHUFFLED=$(wc -l < "$SHUFFLED")
    echo "[gnomad_int] Shuffled (background):  $SHUFFLED_OVERLAP / $N_SHUFFLED"
fi

# -------------------------------------------------------------
# 4. SVtype breakdown of overlapping variants
# -------------------------------------------------------------
echo "[gnomad_int] SV type breakdown of overlaps:"
awk '{print $10}' "$OUT_OVERLAP" | sort | uniq -c | sort -rn | head -10 || true

# -------------------------------------------------------------
# 5. Produce enrichment statistics table
# -------------------------------------------------------------
ENRICH_TABLE="$RESULTS/intersections/enrichment_table_${TIMESTAMP}.tsv"
python3 - <<PYEOF > "$ENRICH_TABLE"
import math

n_delta = $N_DELTA
n_overlap = $N_OVERLAP
n_hifst_overlap = $N_HIFST_OVERLAP
n_gnomad = $(wc -l < "$GNOMAD_BED" 2>/dev/null || echo 0)
n_gnomad_hifst = $(wc -l < "$GNOMAD_HIFST" 2>/dev/null || echo 0)

# Approximate genome size for expected rate
genome_size_bp = 3_000_000_000
avg_delta_len_bp = 300   # typical SV length
expected_rate = (n_gnomad * avg_delta_len_bp) / genome_size_bp
expected_rate_hifst = (n_gnomad_hifst * avg_delta_len_bp) / genome_size_bp

obs_rate = n_overlap / n_delta if n_delta > 0 else 0
obs_rate_hifst = n_hifst_overlap / n_delta if n_delta > 0 else 0

enrichment = obs_rate / expected_rate if expected_rate > 0 else float('inf')
enrichment_hifst = obs_rate_hifst / expected_rate_hifst if expected_rate_hifst > 0 else float('inf')

print("metric\tvalue")
print(f"n_delta_alleles\t{n_delta}")
print(f"n_gnomad_svs\t{n_gnomad}")
print(f"n_gnomad_hifst_svs\t{n_gnomad_hifst}")
print(f"n_overlap_any\t{n_overlap}")
print(f"n_overlap_hifst\t{n_hifst_overlap}")
print(f"overlap_rate_any\t{obs_rate:.4f}")
print(f"overlap_rate_hifst\t{obs_rate_hifst:.4f}")
print(f"expected_rate_any\t{expected_rate:.4f}")
print(f"expected_rate_hifst\t{expected_rate_hifst:.4f}")
print(f"enrichment_any\t{enrichment:.2f}x")
print(f"enrichment_hifst\t{enrichment_hifst:.2f}x")
PYEOF

echo "[gnomad_int] Enrichment table: $ENRICH_TABLE"
cat "$ENRICH_TABLE"
echo "[gnomad_int] Done."
