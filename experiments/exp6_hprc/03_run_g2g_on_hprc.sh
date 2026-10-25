#!/usr/bin/env bash
# ============================================================
# 03_run_g2g_on_hprc.sh - Run G2G on AFR vs EAS pangenomes
#
# Aligns the AFR population pangenome against the EAS pangenome
# using the production G2G binary. Produces:
#   - afr_vs_eas.gaf   alignment in GAF-like format
#   - afr_vs_eas.vcf   delta alleles (population-specific SVs)
#
# Usage:
#   conda activate g2g_eval
#   bash experiments/exp6_hprc/03_run_g2g_on_hprc.sh
# ============================================================
set -euo pipefail

DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
HPRC_DIR="$DATA_ROOT/hprc"
PANGENOME_DIR="$HPRC_DIR/pangenomes"
RESULTS="$DATA_ROOT/results/exp6"
G2G_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
G2G_BIN="$G2G_ROOT/build/release/g2g"
THREADS="${THREADS:-$(nproc)}"
TIMESTAMP=$(date +%Y%m%dT%H%M%S)

mkdir -p "$RESULTS"

echo "[g2g_hprc] Results: $RESULTS"
echo "[g2g_hprc] Threads: $THREADS"

if [[ ! -f "$G2G_BIN" ]]; then
    echo "[g2g_hprc] ERROR: G2G binary not found: $G2G_BIN"
    echo "[g2g_hprc] Run: bash experiments/01_build_all.sh"
    exit 1
fi

# -------------------------------------------------------------
# Find GFA files from pangenome build
# -------------------------------------------------------------
AFR_GFA="$PANGENOME_DIR/AFR/afr_pangenome.gfa"
EAS_GFA="$PANGENOME_DIR/EAS/eas_pangenome.gfa"

# Fallback: look for any .gfa file in the directory
if [[ ! -f "$AFR_GFA" ]]; then
    AFR_GFA=$(find "$PANGENOME_DIR/AFR" -name "*.gfa" | head -1)
fi
if [[ ! -f "$EAS_GFA" ]]; then
    EAS_GFA=$(find "$PANGENOME_DIR/EAS" -name "*.gfa" | head -1)
fi

if [[ -z "$AFR_GFA" ]] || [[ ! -f "$AFR_GFA" ]]; then
    echo "[g2g_hprc] ERROR: AFR GFA not found. Run 02_build_pangenomes.sh first."
    exit 1
fi
if [[ -z "$EAS_GFA" ]] || [[ ! -f "$EAS_GFA" ]]; then
    echo "[g2g_hprc] ERROR: EAS GFA not found. Run 02_build_pangenomes.sh first."
    exit 1
fi

echo "[g2g_hprc] AFR GFA: $AFR_GFA ($(du -sh "$AFR_GFA" | cut -f1))"
echo "[g2g_hprc] EAS GFA: $EAS_GFA ($(du -sh "$EAS_GFA" | cut -f1))"

# -------------------------------------------------------------
# Run G2G
# -------------------------------------------------------------
OUT_GAF="$RESULTS/afr_vs_eas_${TIMESTAMP}.gaf"
OUT_VCF="$RESULTS/afr_vs_eas_${TIMESTAMP}.vcf"
TIME_FILE="$RESULTS/g2g_hprc_time_${TIMESTAMP}.txt"

echo "[g2g_hprc] Running G2G ($THREADS threads)..."
START_EPOCH=$(date +%s)

/usr/bin/time -v -o "$TIME_FILE" \
    "$G2G_BIN" "$AFR_GFA" "$EAS_GFA" \
    --threads "$THREADS" \
    --gaf "$OUT_GAF" \
    --vcf "$OUT_VCF" \
    2>&1 | tee "$RESULTS/g2g_hprc_log_${TIMESTAMP}.txt"

END_EPOCH=$(date +%s)
ELAPSED=$(( END_EPOCH - START_EPOCH ))
PEAK_KB=$(grep "Maximum resident" "$TIME_FILE" 2>/dev/null | grep -o '[0-9]*' | tail -1 || echo 0)

echo "[g2g_hprc] Completed in ${ELAPSED}s  peak_mem=$((PEAK_KB/1024))MB"
echo "g2g_hprc start=$START_EPOCH end=$END_EPOCH elapsed_s=$ELAPSED peak_kb=$PEAK_KB" \
    >> "$RESULTS/timing_${TIMESTAMP}.txt"

# -------------------------------------------------------------
# Summarize output
# -------------------------------------------------------------
if [[ -f "$OUT_GAF" ]]; then
    SNARL_PAIRS=$(grep -c "^" "$OUT_GAF" 2>/dev/null || echo 0)
    echo "[g2g_hprc] Snarl pairs aligned:  $SNARL_PAIRS"
fi

if [[ -f "$OUT_VCF" ]]; then
    DELTA_ALLELES=$(grep -v "^#" "$OUT_VCF" | wc -l)
    echo "[g2g_hprc] Delta alleles (VCF):  $DELTA_ALLELES"
fi

# -------------------------------------------------------------
# Convert VCF to BED for downstream intersection
# -------------------------------------------------------------
OUT_BED="$RESULTS/afr_vs_eas_delta_${TIMESTAMP}.bed"
echo "[g2g_hprc] Converting VCF to BED for bedtools intersection..."

if [[ -f "$OUT_VCF" ]]; then
    python3 - "$OUT_VCF" "$OUT_BED" <<'PYEOF'
import sys, re

in_vcf, out_bed = sys.argv[1], sys.argv[2]
written = 0
with open(in_vcf) as fin, open(out_bed, 'w') as fout:
    for line in fin:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) < 8: continue
        chrom, pos, vid, ref, alt = parts[0], int(parts[1])-1, parts[2], parts[3], parts[4]
        info = dict(kv.split('=',1) for kv in parts[7].split(';') if '=' in kv)
        svtype = info.get('SVTYPE', 'UNK')
        seq    = info.get('SEQ', alt)
        svlen  = len(seq) - len(ref)
        end    = pos + max(1, len(ref))
        # BED: chrom start end id svtype svlen seq
        fout.write(f"{chrom}\t{pos}\t{end}\t{vid}\t{svtype}\t{svlen}\t{seq[:100]}\n")
        written += 1

print(f"[vcf2bed] Wrote {written} records to {out_bed}", file=sys.stderr)
PYEOF
    sort -k1,1 -k2,2n "$OUT_BED" > "${OUT_BED%.bed}_sorted.bed"
    bgzip -@ 4 "${OUT_BED%.bed}_sorted.bed"
    tabix -p bed "${OUT_BED%.bed}_sorted.bed.gz"
    echo "[g2g_hprc] BED (sorted+indexed): ${OUT_BED%.bed}_sorted.bed.gz"
fi

echo "[g2g_hprc] Done. Output files:"
ls -lh "$RESULTS/"*"${TIMESTAMP}"* 2>/dev/null || true
