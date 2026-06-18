#!/usr/bin/env bash
# ============================================================
# 02_build_pangenomes.sh — Build AFR and EAS pangenomes
#
# Uses minigraph-cactus (the HPRC recommended pipeline) to
# build two population-specific pangenome graphs from the
# HPRC assemblies downloaded in download_hprc_assemblies.sh.
#
# Expected runtime: ~48 CPU-hours per pangenome.
# Memory: ~100 GB per pangenome.
#
# Usage:
#   conda activate g2g_eval
#   bash experiments/exp6_hprc/02_build_pangenomes.sh
#
# Set G2G_DATA_ROOT to your data directory (default: $HOME/g2g_data)
# ============================================================
set -euo pipefail

DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
HPRC_DIR="$DATA_ROOT/hprc"
ASSEMBLY_DIR="$HPRC_DIR/assemblies"
METADATA_DIR="$HPRC_DIR/metadata"
PANGENOME_DIR="$HPRC_DIR/pangenomes"
NCPU="${THREADS:-$(nproc)}"
TIMESTAMP=$(date +%Y%m%dT%H%M%S)

mkdir -p "$PANGENOME_DIR"/{AFR,EAS}

echo "[build_pg] Data root:    $HPRC_DIR"
echo "[build_pg] CPUs:         $NCPU"
echo "[build_pg] Timestamp:    $TIMESTAMP"

# Verify cactus is installed
if ! command -v cactus-pangenome &>/dev/null; then
    echo "[build_pg] ERROR: cactus-pangenome not found."
    echo "[build_pg] Install with: conda install -c bioconda cactus"
    exit 1
fi
echo "[build_pg] cactus: $(cactus-pangenome --version 2>&1 | head -1)"

# ─────────────────────────────────────────────────────────────
# Build one pangenome
# ─────────────────────────────────────────────────────────────
build_pangenome() {
    local pop="$1"
    local samples_file="$METADATA_DIR/${pop,,}_samples.txt"
    local out_dir="$PANGENOME_DIR/$pop"
    local seqfile="$out_dir/seqfile.txt"
    local job_store="$out_dir/jobstore"
    local log_file="$out_dir/build_${TIMESTAMP}.log"

    echo ""
    echo "[build_pg] ══ Building $pop pangenome ══"

    if [[ ! -f "$samples_file" ]]; then
        echo "[build_pg] Sample file not found: $samples_file"
        echo "[build_pg] Run 01_stratify_hprc_samples.py first."
        return 1
    fi

    mapfile -t SAMPLES < "$samples_file"
    echo "[build_pg] Samples: ${#SAMPLES[@]}"

    # ── Create cactus-pangenome seqfile (TSV: name  path) ────
    # Format: one line per haplotype
    # <sample>_hap1  /path/to/sample_hap1.fa.gz
    > "$seqfile"

    # Add GRCh38 as the reference (required by cactus-pangenome)
    GRCH38="$DATA_ROOT/reference/GRCh38.fa.gz"
    if [[ -f "$GRCH38" ]]; then
        echo -e "GRCh38\t$GRCH38" >> "$seqfile"
    else
        echo "[build_pg] WARNING: GRCh38 reference not found at $GRCH38"
        echo "[build_pg] Download with: bash experiments/download/download_grch38.sh"
    fi

    for sample in "${SAMPLES[@]}"; do
        hap1="$ASSEMBLY_DIR/$pop/$sample/${sample}_hap1.fa.gz"
        hap2="$ASSEMBLY_DIR/$pop/$sample/${sample}_hap2.fa.gz"
        if [[ -f "$hap1" ]]; then
            echo -e "${sample}_hap1\t$hap1" >> "$seqfile"
        else
            echo "[build_pg] WARNING: missing $hap1"
        fi
        if [[ -f "$hap2" ]]; then
            echo -e "${sample}_hap2\t$hap2" >> "$seqfile"
        else
            echo "[build_pg] WARNING: missing $hap2"
        fi
    done

    echo "[build_pg] Seqfile: $seqfile ($(wc -l < "$seqfile") entries)"

    # ── Run cactus-pangenome ──────────────────────────────────
    START_EPOCH=$(date +%s)
    echo "[build_pg] Starting cactus-pangenome at $(date)..."

    cactus-pangenome \
        "$job_store" \
        "$seqfile" \
        --outDir "$out_dir" \
        --outName "${pop,,}_pangenome" \
        --reference GRCh38 \
        --vcf \
        --gfa \
        --gbz \
        --maxCores "$NCPU" \
        --maxMemory "${SLURM_MEM_PER_NODE:-400}G" \
        --logFile "$log_file" \
        2>&1 | tee -a "$log_file"

    END_EPOCH=$(date +%s)
    ELAPSED=$(( END_EPOCH - START_EPOCH ))
    echo "[build_pg] $pop pangenome complete in ${ELAPSED}s ($(date))"

    # Record timing
    echo "pop=$pop start=$START_EPOCH end=$END_EPOCH elapsed_s=$ELAPSED" \
        >> "$PANGENOME_DIR/build_timing.txt"

    # Verify GFA was produced
    GFA_OUT="$out_dir/${pop,,}_pangenome.gfa"
    if [[ -f "$GFA_OUT" ]]; then
        echo "[build_pg] GFA: $GFA_OUT ($(du -sh "$GFA_OUT" | cut -f1))"
    else
        echo "[build_pg] WARNING: GFA not found at $GFA_OUT"
        ls "$out_dir"/ || true
    fi
}

# ─────────────────────────────────────────────────────────────
# Build both pangenomes (sequential to avoid memory contention)
# ─────────────────────────────────────────────────────────────
build_pangenome AFR
build_pangenome EAS

echo ""
echo "[build_pg] Both pangenomes complete."
echo "[build_pg] Output directory: $PANGENOME_DIR"
echo "[build_pg] GFA files:"
find "$PANGENOME_DIR" -name "*.gfa" -exec ls -lh {} \;
