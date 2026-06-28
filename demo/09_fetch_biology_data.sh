#!/usr/bin/env bash
# 09_fetch_biology_data.sh — stage the public inputs for the REAL Q4 run
# (demo/09_hprc_biology.sh with no --synthetic flag).
#
# Everything here is GRCh37/b37 so coordinates line up end to end:
#   * reference        human_g1k_v37.fasta            (contigs "1".."22", no "chr")
#   * population panel  1000G phase-3 sample->pop map + per-chromosome genotype VCFs
#   * F_ST targets      gnomAD-SV v2.1 sites (Collins 2020 — the paper's citation)
#   * pathogenic sites  ClinVar GRCh37
# Do NOT mix in a "chr"-prefixed GRCh38 file: the contig names must match.
#
# Downloads are large (the reference alone is ~850 MB gz / ~3 GB). Start with one
# small chromosome and widen once the pipeline runs end to end:
#     CHROMS="22"  bash demo/09_fetch_biology_data.sh
#     CHROMS="20 21 22" bash demo/09_fetch_biology_data.sh
#
# This script downloads only; it does not need vg/bcftools (those are used by
# 09_hprc_biology.sh). It cannot run in the dev sandbox (no outbound network there)
# — run it on the cluster. Re-running is safe: existing non-empty files are kept.
set -euo pipefail
DATA="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/data"; mkdir -p "$DATA"; cd "$DATA"
CHROMS="${CHROMS:-22}"
G1K="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp"

get() {  # url [outfile]   — skip if already present and non-empty
    local url="$1" out="${2:-$(basename "$1")}"
    if [ -s "$out" ]; then echo "  have   $out"; return; fi
    echo "  fetch  $out"
    wget -q --show-progress -O "$out" "$url" || { echo "ERROR: download failed: $url" >&2; rm -f "$out"; exit 1; }
}
get_alt() {  # outfile url1 [url2 ...]   — try each mirror until one works
    local out="$1"; shift
    if [ -s "$out" ]; then echo "  have   $out"; return; fi
    local url
    for url in "$@"; do
        echo "  fetch  $out  <- $url"
        if wget -q --show-progress -O "$out" "$url"; then return; fi
        echo "  (mirror failed, trying next)"; rm -f "$out"
    done
    echo "ERROR: all mirrors failed for $out" >&2; exit 1
}
nonempty() { [ -s "$1" ] || { echo "ERROR: $1 is empty or missing" >&2; exit 1; }; }

echo "[1/5] reference (human_g1k_v37) + index + chrom sizes"
get "$G1K/technical/reference/human_g1k_v37.fasta.gz"
[ -s human_g1k_v37.fasta ] || { echo "  gunzip human_g1k_v37.fasta.gz"; gunzip -k human_g1k_v37.fasta.gz; }
if [ ! -s human_g1k_v37.fasta.fai ]; then
    if command -v samtools >/dev/null; then samtools faidx human_g1k_v37.fasta
    else echo "  WARN: samtools not found — vg construct needs human_g1k_v37.fasta.fai;"
         echo "        install it before running 09:  conda install -c bioconda samtools"; fi
fi
if [ -s human_g1k_v37.fasta.fai ]; then
    cut -f1,2 human_g1k_v37.fasta.fai > b37.chrom.sizes
else
    # No samtools yet: write the canonical GRCh37 primary-contig lengths so the
    # enrichment step (which only needs chrom->length) is unblocked.
    printf '%s\t%s\n' 1 249250621 2 243199373 3 198022430 4 191154276 5 180915260 \
        6 171115067 7 159138663 8 146364022 9 141213431 10 135534747 11 135006516 \
        12 133851895 13 115169878 14 107349540 15 102531392 16 90354753 17 81195210 \
        18 78077248 19 59128983 20 63025520 21 48129895 22 51304566 \
        X 155270560 Y 59373566 MT 16569 > b37.chrom.sizes
    echo "  wrote b37.chrom.sizes from canonical GRCh37 lengths"
fi

echo "[2/5] 1000G phase-3 sample -> population panel; derive AFR / EAS lists"
get "$G1K/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" panel.tsv
# columns: sample  pop  super_pop  gender   (super_pop in AFR/AMR/EAS/EUR/SAS)
awk '$3=="AFR"{print $1}' panel.tsv > afr.samples
awk '$3=="EAS"{print $1}' panel.tsv > eas.samples
nonempty afr.samples; nonempty eas.samples
echo "        AFR=$(wc -l < afr.samples)  EAS=$(wc -l < eas.samples) samples"

echo "[3/5] 1000G phase-3 per-chromosome genotype VCFs (CHROMS=$CHROMS)"
for c in $CHROMS; do
    get "$G1K/release/20130502/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \
        "ALL.chr${c}.phase3.genotypes.vcf.gz"
done

echo "[4/5] gnomAD-SV v2.1 sites (GRCh37) — GCS, with AWS S3 mirror as fallback"
get_alt gnomad_v2.1_sv.sites.vcf.gz \
    "https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz" \
    "https://gnomad-public-us-east-1.s3.amazonaws.com/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz"

echo "[5/5] ClinVar (GRCh37)"
get "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz" clinvar_GRCh37.vcf.gz

echo
echo "gnomAD-SV per-population AF fields — confirm the names for --afr-field/--eas-field:"
zcat gnomad_v2.1_sv.sites.vcf.gz 2>/dev/null | grep -m 200 '^##INFO' | grep -iE 'id=(afr|eas)' \
    || echo "  (none matched 'afr/eas' — inspect: zcat gnomad_v2.1_sv.sites.vcf.gz | grep '^##INFO')"

cat <<EOF

Done. Inputs are in $DATA. Now run the real Q4 pipeline (start with CHROMS="$CHROMS"):

  REF_FA=$DATA/human_g1k_v37.fasta \\
  PANEL_TMPL=$DATA/ALL.chr{}.phase3.genotypes.vcf.gz \\
  SAMPLES_A=$DATA/afr.samples  SAMPLES_B=$DATA/eas.samples \\
  GNOMAD_SV_VCF=$DATA/gnomad_v2.1_sv.sites.vcf.gz \\
  CLINVAR_VCF=$DATA/clinvar_GRCh37.vcf.gz \\
  CHROM_SIZES=$DATA/b37.chrom.sizes \\
  CHROMS="$CHROMS"  AFR_FIELD=AFR_AF  EAS_FIELD=EAS_AF \\
  bash demo/09_hprc_biology.sh

(If the header above shows different AF field names, set AFR_FIELD/EAS_FIELD to match.)
EOF
