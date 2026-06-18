#!/usr/bin/env bash
# ============================================================
# 05_intersect_clinvar.sh — Intersect G2G shared alleles with ClinVar
#
# Tests whether the *shared* G2G alleles (present in both AFR
# and EAS pangenomes) contain known pathogenic variants from
# ClinVar — demonstrating cross-population conservation of
# clinically relevant SVs.
#
# Hypothesis:
#   Shared alleles should be enriched for pathogenic/LP ClinVar
#   SVs (variants conserved because they affect essential biology).
#
# Usage:
#   bash experiments/exp6_hprc/05_intersect_clinvar.sh
# ============================================================
set -euo pipefail

DATA_ROOT="${G2G_DATA_ROOT:-/data/fat/g2g}"
RESULTS="$DATA_ROOT/results/exp6"
CLINVAR_DIR="$DATA_ROOT/clinvar"
TIMESTAMP=$(date +%Y%m%dT%H%M%S)

mkdir -p "$RESULTS/intersections"

# ─────────────────────────────────────────────────────────────
# Find G2G alignment output (GAF file has the matched alleles)
# ─────────────────────────────────────────────────────────────
GAF_FILE=$(ls -t "$RESULTS"/afr_vs_eas_*.gaf 2>/dev/null | head -1)
if [[ -z "$GAF_FILE" ]]; then
    echo "[clinvar_int] ERROR: No GAF file found in $RESULTS"
    exit 1
fi
echo "[clinvar_int] GAF: $GAF_FILE"

CLINVAR_PATH="$CLINVAR_DIR/clinvar_pathogenic_svs.bed"
CLINVAR_ALL="$CLINVAR_DIR/clinvar_all_svs.bed"

if [[ ! -f "$CLINVAR_PATH" ]]; then
    echo "[clinvar_int] ERROR: ClinVar BED not found: $CLINVAR_PATH"
    exit 1
fi

# ─────────────────────────────────────────────────────────────
# 1. Extract shared allele coordinates from GAF
#    GAF format: each line = one snarl alignment
#    Shared alleles are the matched allele_matches (not delta)
# ─────────────────────────────────────────────────────────────
SHARED_BED="$RESULTS/intersections/shared_alleles_${TIMESTAMP}.bed"
echo "[clinvar_int] Extracting shared allele coordinates from GAF..."

python3 - "$GAF_FILE" "$SHARED_BED" <<'PYEOF'
import sys, re

gaf_file, out_bed = sys.argv[1], sys.argv[2]
written = 0

with open(gaf_file) as fin, open(out_bed, 'w') as fout:
    for line in fin:
        line = line.strip()
        if not line or line.startswith('#'): continue
        # GAF-like format from g2g:
        # Each alignment line has tab-separated fields.
        # We look for coordinate annotations in INFO-like fields.
        parts = line.split('\t')
        # Field parsing depends on exact GAF writer format in io.cpp
        # Extract position from path field (first token = snarl_id "chr:pos-end")
        if len(parts) < 5: continue
        snarl_id = parts[0]  # format: "chrom:src_pos>snk_pos" or similar
        m = re.match(r'^(\w+):(\d+)>(\d+)', snarl_id)
        if m:
            chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
            score = parts[4] if len(parts) > 4 else "0"
            fout.write(f"{chrom}\t{start}\t{end}\t{snarl_id}\tshared\t{score}\n")
            written += 1

print(f"[gaf2bed] Wrote {written} shared allele regions to {out_bed}", file=sys.stderr)
PYEOF

if [[ ! -s "$SHARED_BED" ]]; then
    echo "[clinvar_int] WARNING: shared BED is empty — GAF format may not include coordinates."
    echo "[clinvar_int] Attempting fallback: extract from VCF complement..."
    # Alternative: use the delta VCF to define shared regions as genome-wide minus delta
    DELTA_VCF=$(ls -t "$RESULTS"/afr_vs_eas_*.vcf 2>/dev/null | head -1)
    if [[ -n "$DELTA_VCF" ]]; then
        python3 - "$DELTA_VCF" "$SHARED_BED" <<'PYEOF2'
import sys
in_vcf, out_bed = sys.argv[1], sys.argv[2]
# A shared allele is one that appears in the matched section of the VCF
# (SHARED=1 annotation)
written = 0
with open(in_vcf) as fin, open(out_bed, 'w') as fout:
    for line in fin:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) < 8: continue
        info = dict(kv.split('=',1) for kv in parts[7].split(';') if '=' in kv)
        if info.get('SHARED','0') == '1':
            chrom, pos = parts[0], int(parts[1])-1
            ref = parts[3]
            end = pos + max(1, len(ref))
            fout.write(f"{chrom}\t{pos}\t{end}\t{parts[2]}\tshared\t{info.get('SCORE','0')}\n")
            written += 1
print(f"[vcf2bed] Wrote {written} shared records", file=sys.stderr)
PYEOF2
    fi
fi

# Sort and index
sort -k1,1 -k2,2n "$SHARED_BED" > "${SHARED_BED%.bed}_sorted.bed"
bgzip -@ 4 "${SHARED_BED%.bed}_sorted.bed"
tabix -p bed "${SHARED_BED%.bed}_sorted.bed.gz"
SHARED_BED_GZ="${SHARED_BED%.bed}_sorted.bed.gz"
echo "[clinvar_int] Shared BED: $SHARED_BED_GZ ($(zcat "$SHARED_BED_GZ" | wc -l) regions)"

# ─────────────────────────────────────────────────────────────
# 2. Intersect shared alleles with ClinVar pathogenic SVs
# ─────────────────────────────────────────────────────────────
OUT_PATHOGENIC="$RESULTS/intersections/shared_vs_clinvar_pathogenic_${TIMESTAMP}.bed"
echo "[clinvar_int] Intersecting with ClinVar pathogenic SVs..."

bedtools intersect \
    -a "$SHARED_BED_GZ" \
    -b "$CLINVAR_PATH" \
    -wa -wb \
    > "$OUT_PATHOGENIC"

N_SHARED=$(zcat "$SHARED_BED_GZ" | wc -l)
N_PATH_OVERLAP=$(wc -l < "$OUT_PATHOGENIC")
echo "[clinvar_int] Shared alleles:         $N_SHARED"
echo "[clinvar_int] ClinVar path. overlaps: $N_PATH_OVERLAP"

# ─────────────────────────────────────────────────────────────
# 3. Intersect with all ClinVar SVs (broader)
# ─────────────────────────────────────────────────────────────
OUT_ALL_CLINVAR="$RESULTS/intersections/shared_vs_clinvar_all_${TIMESTAMP}.bed"
bedtools intersect \
    -a "$SHARED_BED_GZ" \
    -b "$CLINVAR_ALL" \
    -wa -wb \
    > "$OUT_ALL_CLINVAR"
N_ALL_CLINVAR_OVERLAP=$(wc -l < "$OUT_ALL_CLINVAR")
echo "[clinvar_int] ClinVar all SVs overlaps: $N_ALL_CLINVAR_OVERLAP"

# ─────────────────────────────────────────────────────────────
# 4. Top overlapping ClinVar conditions
# ─────────────────────────────────────────────────────────────
if [[ -s "$OUT_PATHOGENIC" ]]; then
    echo ""
    echo "[clinvar_int] Top ClinVar conditions overlapping shared alleles:"
    awk '{print $NF}' "$OUT_PATHOGENIC" | sort | uniq -c | sort -rn | head -15
fi

# ─────────────────────────────────────────────────────────────
# 5. Summary statistics table
# ─────────────────────────────────────────────────────────────
CLINVAR_TABLE="$RESULTS/intersections/clinvar_summary_${TIMESTAMP}.tsv"
python3 - <<PYEOF > "$CLINVAR_TABLE"
n_shared = $N_SHARED
n_pathogenic = $N_PATH_OVERLAP
n_all_clinvar = $N_ALL_CLINVAR_OVERLAP
n_clinvar_pathogenic_total = $(wc -l < "$CLINVAR_PATH" 2>/dev/null || echo 1)
n_clinvar_all_total = $(wc -l < "$CLINVAR_ALL" 2>/dev/null || echo 1)

pathogenic_recall = n_pathogenic / n_clinvar_pathogenic_total if n_clinvar_pathogenic_total > 0 else 0
any_recall = n_all_clinvar / n_clinvar_all_total if n_clinvar_all_total > 0 else 0
pathogenic_rate = n_pathogenic / n_shared if n_shared > 0 else 0

print("metric\tvalue")
print(f"n_shared_alleles\t{n_shared}")
print(f"n_clinvar_pathogenic_svs\t{n_clinvar_pathogenic_total}")
print(f"n_clinvar_all_svs\t{n_clinvar_all_total}")
print(f"n_shared_overlapping_pathogenic\t{n_pathogenic}")
print(f"n_shared_overlapping_any_clinvar\t{n_all_clinvar}")
print(f"pathogenic_recall\t{pathogenic_recall:.4f}")
print(f"any_clinvar_recall\t{any_recall:.4f}")
print(f"pathogenic_overlap_rate\t{pathogenic_rate:.4f}")
PYEOF

echo "[clinvar_int] ClinVar summary table:"
cat "$CLINVAR_TABLE"
echo "[clinvar_int] Done."
