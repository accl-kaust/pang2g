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
# 3. Resolve real per-haplotype fasta URLs from the HPRC index
# ─────────────────────────────────────────────────────────────
# The assembly fasta locations are NOT a guessable S3 path — they are listed
# (as s3:// URLs) inside the index file we just downloaded, in columns whose
# names contain "hap1"/"hap2" and "fasta". We parse them into a simple
# sample -> hap1_url, hap2_url map and convert s3:// to https.
INDEX="$METADATA_DIR/hprc_v2_index.tsv"
URLMAP="$METADATA_DIR/assembly_urls.tsv"

python3 - "$INDEX" "$URLMAP" <<'PYEOF'
import sys, csv
idx, out = sys.argv[1], sys.argv[2]
with open(idx, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    cols = reader.fieldnames or []

    def find(hap):
        # prefer a column that is clearly an AWS/S3 fasta URL for this hap
        for want in (lambda c: hap in c and 'fasta' in c and ('aws' in c or 's3' in c),
                     lambda c: hap in c and 'fasta' in c):
            for c in cols:
                if want(c.lower()):
                    return c
        return None

    c1, c2 = find('hap1'), find('hap2')
    scol = next((c for c in cols if c.lower() in
                 ('sample', 'sample_id', 'sampleid')), cols[0] if cols else None)

    sys.stderr.write(f"[hprc] index columns: {cols}\n")
    sys.stderr.write(f"[hprc] using sample='{scol}' hap1='{c1}' hap2='{c2}'\n")
    if not (c1 and c2 and scol):
        sys.stderr.write("[hprc] ERROR: could not locate fasta-URL columns in index\n")
        sys.exit(2)

    n = 0
    with open(out, 'w') as o:
        for row in reader:
            s = (row.get(scol) or '').strip()
            u1 = (row.get(c1) or '').strip()
            u2 = (row.get(c2) or '').strip()
            if s:
                o.write(f"{s}\t{u1}\t{u2}\n"); n += 1
    sys.stderr.write(f"[hprc] wrote {n} sample URL rows -> {out}\n")
PYEOF

# ─────────────────────────────────────────────────────────────
# 4. Download assemblies (both haplotypes per sample)
# ─────────────────────────────────────────────────────────────
to_https() {
    local u="$1"
    if [[ "$u" == s3://* ]]; then
        local rest="${u#s3://}"
        echo "https://${rest%%/*}.s3.amazonaws.com/${rest#*/}"
    else
        echo "$u"
    fi
}
export -f to_https

download_sample() {
    local sample="$1"
    local pop="$2"
    local out_dir="$ASSEMBLY_DIR/$pop/$sample"
    mkdir -p "$out_dir"

    local row; row=$(grep -m1 -P "^${sample}\t" "$URLMAP" || true)
    if [[ -z "$row" ]]; then
        echo "[hprc] WARNING: $sample not found in index URL map"
        return
    fi
    local u1 u2
    u1=$(cut -f2 <<< "$row"); u2=$(cut -f3 <<< "$row")

    for hap in 1 2; do
        local url; [[ "$hap" == "1" ]] && url=$(to_https "$u1") || url=$(to_https "$u2")
        local out_file="$out_dir/${sample}_hap${hap}.fa.gz"

        if [[ -f "$out_file" ]] && [[ $(stat -c%s "$out_file") -gt 1000000 ]]; then
            echo "[hprc] SKIP (exists): $out_file"
            continue
        fi
        if [[ -z "$url" ]]; then
            echo "[hprc] WARNING: no hap$hap URL for $sample"
            continue
        fi

        echo "[hprc] Downloading: $sample hap$hap ($pop)"
        wget -q --show-progress --retry-connrefused --tries=5 \
             --waitretry=30 -O "$out_file" "$url" || \
            echo "[hprc] WARNING: could not download $sample hap$hap ($url)"
    done
}

export -f download_sample
export ASSEMBLY_DIR URLMAP

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
# 5. Verify downloads
# ─────────────────────────────────────────────────────────────
echo ""
echo "[hprc] Verifying downloads..."
AFR_COUNT=$(find "$ASSEMBLY_DIR/AFR" -name "*.fa.gz" -size +1M | wc -l)
EAS_COUNT=$(find "$ASSEMBLY_DIR/EAS" -name "*.fa.gz" -size +1M | wc -l)
echo "[hprc] AFR assemblies found: $AFR_COUNT (expected: $((N_PER_POP * 2)))"
echo "[hprc] EAS assemblies found: $EAS_COUNT (expected: $((N_PER_POP * 2)))"
echo "[hprc] Total disk usage: $(du -sh "$ASSEMBLY_DIR")"

if [[ "$AFR_COUNT" -eq 0 && "$EAS_COUNT" -eq 0 ]]; then
    echo "[hprc] ERROR: no assemblies downloaded. Inspect the index header" >&2
    echo "[hprc]        and URL map: $URLMAP" >&2
    echo "[hprc]        First rows of the index:" >&2
    head -2 "$INDEX" >&2 || true
    exit 1
fi
