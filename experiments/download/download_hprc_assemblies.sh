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

# Download the authoritative 1000 Genomes super-population panel so we can
# classify each HPRC sample by ancestry instead of trusting a hardcoded list.
# Try several mirrors; validate the result actually looks like the panel
# (a real sample row mentioning a super-population) rather than an empty file
# or an HTML error page — otherwise classification silently yields 0 samples.
PANEL="$METADATA_DIR/1kgp_pops.panel"
PANEL_URLS=(
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
    "https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
    "https://s3.amazonaws.com/1000genomes/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
)
echo "[hprc] Downloading 1KGP super-population panel..."
PANEL_OK=0
for u in "${PANEL_URLS[@]}"; do
    if wget -q -O "$PANEL" "$u" 2>/dev/null || curl -fsSL -o "$PANEL" "$u" 2>/dev/null; then
        if grep -qE '(AFR|EAS)' "$PANEL" 2>/dev/null && \
           [[ $(wc -l < "$PANEL") -gt 1000 ]]; then
            PANEL_OK=1
            echo "[hprc] panel OK ($(wc -l < "$PANEL") rows) from $u"
            break
        fi
    fi
done
if [[ "$PANEL_OK" -eq 0 ]]; then
    echo "[hprc] ERROR: could not fetch a valid 1KGP panel from any mirror." >&2
    echo "[hprc]        Population classification needs it. Check login-node" >&2
    echo "[hprc]        network, or place the panel at: $PANEL" >&2
    exit 1
fi

echo "[hprc] Metadata written to $METADATA_DIR"

# ─────────────────────────────────────────────────────────────
# 2+3. Resolve fasta URLs AND select a BALANCED AFR/EAS sample set
# ─────────────────────────────────────────────────────────────
# The HPRC Year-1 assembly set is heavily skewed away from EAS — only ~6 East
# Asian samples exist in the whole index. Rather than hardcode (fragile, and our
# earlier guesses were wrong: NA18*/NA19* are actually YRI), we join the index
# against the 1KGP panel to get each sample's super-population, then pick an
# EQUAL number from AFR and EAS (capped by EAS availability and N_PER_POP) for a
# size-matched comparison. Fasta URLs live (as s3:// URLs) inside the index in
# columns named like "hap1_*_fasta"; we parse and convert them to https below.
INDEX="$METADATA_DIR/hprc_v2_index.tsv"
URLMAP="$METADATA_DIR/assembly_urls.tsv"

python3 - "$INDEX" "$PANEL" "$URLMAP" "$METADATA_DIR" "$N_PER_POP" <<'PYEOF'
import sys, csv, os
idx, panel, out, metadir, npp = sys.argv[1:6]
npp = int(npp)

# Build sample -> super-population map from the 1KGP panel.
superpop = {}
try:
    with open(panel) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[0] != "sample":
                superpop[p[0]] = p[2]
except FileNotFoundError:
    sys.stderr.write("[hprc] WARNING: 1KGP panel missing; GIAB supplement only\n")

# GIAB / non-1KGP samples the panel does not cover.
superpop.setdefault("HG005", "EAS")   # GIAB Han Chinese son (NA24631)
sys.stderr.write(f"[hprc] super-pop entries loaded: {len(superpop)}\n")

with open(idx, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    cols = reader.fieldnames or []

    def find(hap):
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

    rows = {}
    for row in reader:
        s = (row.get(scol) or '').strip()
        if s:
            rows[s] = ((row.get(c1) or '').strip(), (row.get(c2) or '').strip())

# Full sample -> URL map (all samples in the index).
with open(out, 'w') as o:
    for s, (u1, u2) in rows.items():
        o.write(f"{s}\t{u1}\t{u2}\n")
sys.stderr.write(f"[hprc] wrote {len(rows)} sample URL rows -> {out}\n")

# Classify available samples by super-population, then pick a balanced set.
avail = {"AFR": [], "EAS": []}
for s in rows:
    sp = superpop.get(s)
    if sp in avail:
        avail[sp].append(s)

n = min(len(avail["AFR"]), len(avail["EAS"]), npp)
sys.stderr.write(f"[hprc] available AFR={len(avail['AFR'])} "
                 f"EAS={len(avail['EAS'])} -> balanced n={n} per pop\n")
if n == 0:
    sys.stderr.write("[hprc] ERROR: could not select a balanced AFR/EAS set\n")
    sys.exit(3)

for pop in ("AFR", "EAS"):
    sel = sorted(avail[pop])[:n]
    with open(os.path.join(metadir, f"{pop.lower()}_samples.txt"), 'w') as o:
        o.write("\n".join(sel) + "\n")
    sys.stderr.write(f"[hprc] {pop} selected ({len(sel)}): {sel}\n")
PYEOF

# Read the balanced selections the Python step just wrote.
mapfile -t AFR_SAMPLES < "$METADATA_DIR/afr_samples.txt"
mapfile -t EAS_SAMPLES < "$METADATA_DIR/eas_samples.txt"
echo "[hprc] AFR samples (${#AFR_SAMPLES[@]}): ${AFR_SAMPLES[*]}"
echo "[hprc] EAS samples (${#EAS_SAMPLES[@]}): ${EAS_SAMPLES[*]}"

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

# Download AFR samples in parallel (lists are already balanced + size-capped)
echo "[hprc] Downloading AFR assemblies (parallel, $NCPU jobs)..."
mkdir -p "$ASSEMBLY_DIR/AFR"
printf '%s\n' "${AFR_SAMPLES[@]}" | \
    xargs -P "$NCPU" -I{} bash -c 'download_sample "$@"' _ {} AFR

# Download EAS samples in parallel
echo "[hprc] Downloading EAS assemblies (parallel, $NCPU jobs)..."
mkdir -p "$ASSEMBLY_DIR/EAS"
printf '%s\n' "${EAS_SAMPLES[@]}" | \
    xargs -P "$NCPU" -I{} bash -c 'download_sample "$@"' _ {} EAS

# ─────────────────────────────────────────────────────────────
# 5. Verify downloads
# ─────────────────────────────────────────────────────────────
echo ""
echo "[hprc] Verifying downloads..."
AFR_COUNT=$(find "$ASSEMBLY_DIR/AFR" -name "*.fa.gz" -size +1M | wc -l)
EAS_COUNT=$(find "$ASSEMBLY_DIR/EAS" -name "*.fa.gz" -size +1M | wc -l)
echo "[hprc] AFR assemblies found: $AFR_COUNT (expected: $((${#AFR_SAMPLES[@]} * 2)))"
echo "[hprc] EAS assemblies found: $EAS_COUNT (expected: $((${#EAS_SAMPLES[@]} * 2)))"
echo "[hprc] Total disk usage: $(du -sh "$ASSEMBLY_DIR")"

if [[ "$AFR_COUNT" -eq 0 && "$EAS_COUNT" -eq 0 ]]; then
    echo "[hprc] ERROR: no assemblies downloaded. Inspect the index header" >&2
    echo "[hprc]        and URL map: $URLMAP" >&2
    echo "[hprc]        First rows of the index:" >&2
    head -2 "$INDEX" >&2 || true
    exit 1
fi
