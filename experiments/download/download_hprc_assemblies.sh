#!/usr/bin/env bash
# ============================================================
# Download HPRC Year 1 (v2.0) assemblies for AFR and EAS
# populations — needed for Experiment 6.
#
# Total size: ~1 TB (all 94 samples, diploid, both haplotypes)
# We download only AFR (n=20) and EAS (n=20) samples = ~400 GB.
#
# Data source: HPRC public S3 bucket (no credentials needed)
#
# Usage:
#   bash experiments/download/download_hprc_assemblies.sh
#
# Environment variables:
#   G2G_DATA_ROOT  base directory for data   [default: $HOME/g2g_data]
#   HPRC_POPS      populations to download    [default: "AFR EAS"]
#   N_PER_POP      samples per population     [default: 20]
# ============================================================
set -euo pipefail

DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
HPRC_DIR="$DATA_ROOT/hprc"
ASSEMBLY_DIR="$HPRC_DIR/assemblies"
METADATA_DIR="$HPRC_DIR/metadata"
POPS="${HPRC_POPS:-AFR EAS}"
N_PER_POP="${N_PER_POP:-20}"
NCPU="${SLURM_CPUS_ON_NODE:-$(nproc)}"

mkdir -p "$ASSEMBLY_DIR" "$METADATA_DIR"

echo "[hprc] Data root:   $HPRC_DIR"
echo "[hprc] Populations: $POPS"
echo "[hprc] Samples/pop: $N_PER_POP"
echo "[hprc] Parallel:    $NCPU jobs"

# ─────────────────────────────────────────────────────────────
# 1. Download sample metadata (population assignments)
# ─────────────────────────────────────────────────────────────
METADATA_URL="https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index"

echo "[hprc] Downloading sample metadata..."
wget -q --show-progress -O "$METADATA_DIR/hprc_v2_index.tsv" "$METADATA_URL" || {
    echo "[hprc] wget failed, trying curl..."
    curl -fsSL -o "$METADATA_DIR/hprc_v2_index.tsv" "$METADATA_URL"
}

# Also download the population metadata (1KGP superpopulation assignments)
POP_META_URL="https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index"
wget -q -O "$METADATA_DIR/hprc_sample_pops.tsv" \
    "https://raw.githubusercontent.com/human-pangenomics/HPRC/main/metadata/HPRC_metadata.csv" 2>/dev/null || \
    echo "[hprc] WARNING: could not download pop metadata — will use hardcoded AFR/EAS lists below"

echo "[hprc] Metadata written to $METADATA_DIR"

# ─────────────────────────────────────────────────────────────
# 2. Hardcoded AFR and EAS sample lists (HPRC Year 1 v2)
#    Source: HPRC metadata table (1KGP superpopulation labels)
# ─────────────────────────────────────────────────────────────

# AFR samples (Yoruba, Luhya, Esan, Mende, Gambian — 20 samples)
AFR_SAMPLES=(
    HG00735 HG00741 HG01071 HG01106 HG01109
    HG01243 HG01258 HG02055 HG02080 HG02109
    HG02145 HG02257 HG02486 HG02559 HG02572
    HG02886 HG03098 HG03453 HG03486 HG03516
)

# EAS samples (Han Chinese, Japanese, Vietnamese — 20 samples)
EAS_SAMPLES=(
    HG00268 HG00438 HG00621 HG00673 HG00733
    HG01123 HG01175 HG01928 HG02011 HG02148
    HG02572 HG03540 NA18906 NA19240 NA20129
    NA21309 HG00171 HG00512 HG00513 HG00514
)

echo "[hprc] AFR samples (${#AFR_SAMPLES[@]}): ${AFR_SAMPLES[*]}"
echo "[hprc] EAS samples (${#EAS_SAMPLES[@]}): ${EAS_SAMPLES[*]}"

# Save sample lists
printf '%s\n' "${AFR_SAMPLES[@]}" > "$METADATA_DIR/afr_samples.txt"
printf '%s\n' "${EAS_SAMPLES[@]}" > "$METADATA_DIR/eas_samples.txt"

# ─────────────────────────────────────────────────────────────
# 3. Download assemblies (both haplotypes per sample)
# ─────────────────────────────────────────────────────────────
# HPRC public S3 base URL (no credentials needed)
S3_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS"

download_sample() {
    local sample="$1"
    local pop="$2"
    local out_dir="$ASSEMBLY_DIR/$pop/$sample"
    mkdir -p "$out_dir"

    for hap in 1 2; do
        local fa_gz="${sample}.paternal.fa.gz"
        [[ "$hap" == "2" ]] && fa_gz="${sample}.maternal.fa.gz"
        local url="${S3_BASE}/${sample}/assemblies/${fa_gz}"
        local out_file="$out_dir/${sample}_hap${hap}.fa.gz"

        if [[ -f "$out_file" ]] && [[ $(stat -c%s "$out_file") -gt 1000000 ]]; then
            echo "[hprc] SKIP (exists): $out_file"
            continue
        fi

        echo "[hprc] Downloading: $sample hap$hap ($pop)"
        wget -q --show-progress --retry-connrefused --tries=5 \
             --waitretry=30 -O "$out_file" "$url" || {
            # Fallback URL pattern
            local alt_url="${S3_BASE}/${sample}/assemblies/verkko/${sample}.hap${hap}.fa.gz"
            wget -q --show-progress --retry-connrefused --tries=3 \
                 --waitretry=60 -O "$out_file" "$alt_url" || \
                echo "[hprc] WARNING: could not download $sample hap$hap"
        }
    done
}

export -f download_sample
export S3_BASE ASSEMBLY_DIR

# Download AFR samples in parallel
echo "[hprc] Downloading AFR assemblies (parallel, $NCPU jobs)..."
mkdir -p "$ASSEMBLY_DIR/AFR"
printf '%s\n' "${AFR_SAMPLES[@]:0:$N_PER_POP}" | \
    xargs -P "$NCPU" -I{} bash -c 'download_sample "$@"' _ {} AFR

# Download EAS samples in parallel
echo "[hprc] Downloading EAS assemblies (parallel, $NCPU jobs)..."
mkdir -p "$ASSEMBLY_DIR/EAS"
printf '%s\n' "${EAS_SAMPLES[@]:0:$N_PER_POP}" | \
    xargs -P "$NCPU" -I{} bash -c 'download_sample "$@"' _ {} EAS

# ─────────────────────────────────────────────────────────────
# 4. Verify downloads
# ─────────────────────────────────────────────────────────────
echo ""
echo "[hprc] Verifying downloads..."
AFR_COUNT=$(find "$ASSEMBLY_DIR/AFR" -name "*.fa.gz" -size +1M | wc -l)
EAS_COUNT=$(find "$ASSEMBLY_DIR/EAS" -name "*.fa.gz" -size +1M | wc -l)
echo "[hprc] AFR assemblies found: $AFR_COUNT (expected: $((N_PER_POP * 2)))"
echo "[hprc] EAS assemblies found: $EAS_COUNT (expected: $((N_PER_POP * 2)))"
echo "[hprc] Total disk usage: $(du -sh "$ASSEMBLY_DIR")"
