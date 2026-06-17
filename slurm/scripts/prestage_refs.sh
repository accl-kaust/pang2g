#!/bin/bash
# =============================================================================
# prestage_refs.sh — Download + prepare exp6 reference/truth data (LOGIN NODE).
# =============================================================================
# exp6 needs THREE reference files that the assembly pre-stage does not cover.
# This script downloads them from their canonical public sources and converts
# them into the EXACT BED layouts exp6's Stage 5 expects:
#
#   $REF_DIR/GRCh38.fa            reference genome (cactus anchor + lift coords)
#   $REF_DIR/gnomad_sv_v4.bed     BED: chrom  start  end  id  fst   (fst in col 5)
#   $REF_DIR/clinvar_sv.bed       BED: chrom  start  end  id
#
# SOURCES (verified public, GRCh38-native):
#   GRCh38   : NCBI analysis set (no-alt), same lineage HPRC/1000G use.
#   gnomAD-SV: gnomAD v4.1 SV sites BED (GCS public bucket).
#   ClinVar  : dbVar nstd102 (ClinVar clinical SVs), GRCh38 VCF.
#
# The gnomAD "fst" column: gnomAD-SV v4 does not publish a per-pair F_ST, so we
# derive an AFR-vs-EAS differentiation proxy from the per-population allele
# frequencies (AF_afr, AF_eas) carried in the BED:
#     diff = |AF_afr - AF_eas|        (0..1; higher = more differentiated)
# This is a transparent stand-in for F_ST that matches the AFR/EAS panels the
# experiment compares. Documented here so it is not mistaken for true Weir-
# Cockerham F_ST. Override the column/source if you have a real F_ST track.
#
# RUN (login node, has internet):
#     bash slurm/scripts/prestage_refs.sh
# Idempotent: existing, non-empty outputs are skipped.
#
# CONFIG (env overrides):
#   G2G_REPO, G2G_DATA_ROOT  (as in prestage_hprc.sh)
#   G2G_GRCH38_URL, G2G_GNOMAD_URL, G2G_CLINVAR_URL  (override any source)
# =============================================================================
set -euo pipefail

if [[ -n "${G2G_REPO:-}" ]]; then REPO="$G2G_REPO"
else REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; fi

DATA_ROOT="${G2G_DATA_ROOT:-$REPO/data/hprc}"
REF_DIR="$DATA_ROOT/reference"
DL_DIR="$REF_DIR/_downloads"
mkdir -p "$REF_DIR" "$DL_DIR"

GRCH38_URL="${G2G_GRCH38_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz}"
GNOMAD_URL="${G2G_GNOMAD_URL:-https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.bed.gz}"
CLINVAR_URL="${G2G_CLINVAR_URL:-https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd102.GRCh38.variant_call.vcf.gz}"

GRCH38="$REF_DIR/GRCh38.fa"
GNOMAD_BED="$REF_DIR/gnomad_sv_v4.bed"
CLINVAR_BED="$REF_DIR/clinvar_sv.bed"
MIN_SV_LEN="${G2G_MIN_SV_LEN:-50}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }

command -v wget >/dev/null 2>&1 || { log "FATAL: wget not found."; exit 1; }

# Quick connectivity preflight (fail fast on an air-gapped node).
log "preflight: checking internet ..."
if ! wget --timeout=15 --tries=1 --spider -q "https://ftp.ncbi.nlm.nih.gov/" 2>/dev/null \
   && ! { command -v curl >/dev/null 2>&1 && curl -sS --max-time 15 -I "https://ftp.ncbi.nlm.nih.gov/" >/dev/null 2>&1; }; then
    log "FATAL: no outbound internet from this machine ($(hostname))."
    log "       Run this on a machine that has internet, then rsync \$REF_DIR over:"
    log "         rsync -avP $REF_DIR/  <user>@mcmgt01:$REF_DIR/"
    exit 1
fi
log "preflight OK"

fetch() {  # fetch <url> <dest>
    local url="$1" dest="$2"
    if [[ -s "$dest" ]]; then log "  cached: $(basename "$dest")"; return 0; fi
    log "  downloading $(basename "$dest")"
    wget --continue --tries=3 --waitretry=15 --timeout=120 --show-progress -O "$dest" "$url" \
        || { log "  FAILED: $url"; rm -f "$dest"; return 1; }
}

# ── 1. GRCh38 reference FASTA ────────────────────────────────────────────────
log "=== 1/3: GRCh38 reference ==="
if [[ -s "$GRCH38" ]]; then
    log "  GRCh38.fa already present (skip)"
else
    fetch "$GRCH38_URL" "$DL_DIR/grch38.fna.gz"
    log "  decompressing -> $GRCH38"
    gunzip -c "$DL_DIR/grch38.fna.gz" > "$GRCH38"
    if command -v samtools >/dev/null 2>&1; then
        samtools faidx "$GRCH38" && log "  indexed (.fai created)"
    else
        log "  NOTE: samtools not found — .fai not created (exp6 will create it if needed)"
    fi
fi

# ── 2. gnomAD-SV v4.1 -> BED (chrom start end id fst) ────────────────────────
log "=== 2/3: gnomAD-SV v4.1 ==="
if [[ -s "$GNOMAD_BED" ]]; then
    log "  gnomad_sv_v4.bed already present (skip)"
else
    fetch "$GNOMAD_URL" "$DL_DIR/gnomad_sv.bed.gz"
    log "  converting -> $GNOMAD_BED (fst = |AF_afr - AF_eas|)"
    # The gnomAD v4.1 SV sites BED is a headered TSV. We locate columns by NAME
    # from the header line so this is robust to column re-ordering: chrom,start,
    # end, name, AF_afr, AF_eas. SVLEN (if present) is used for the size floor.
    zcat "$DL_DIR/gnomad_sv.bed.gz" | python3 -c '
import sys
hdr=None
idx={}
def find(cols, *names):
    for n in names:
        if n in cols: return cols.index(n)
    return -1
out=sys.stdout
for line in sys.stdin:
    line=line.rstrip("\n")
    if not line: continue
    if hdr is None:
        # header line: starts with #chrom or chrom
        cols=line.lstrip("#").split("\t")
        hdr=cols
        ic=find(cols,"chrom","#chrom")
        istart=find(cols,"start","pos","chromStart")
        iend=find(cols,"end","chromEnd")
        iname=find(cols,"name","ID","svid")
        iafr=find(cols,"AF_afr","afr_AF","AFR_AF")
        ieas=find(cols,"AF_eas","eas_AF","EAS_AF")
        isvlen=find(cols,"SVLEN","svlen")
        # fall back to positional if header is absent / unrecognised
        if ic<0 or istart<0 or iend<0:
            # treat this first line as data with positional cols
            ic,istart,iend,iname=0,1,2,3; iafr=ieas=isvlen=-1
            hdr=[]  # so we reprocess this line as data
        else:
            continue
    f=line.split("\t")
    try:
        chrom=f[ic]; start=int(float(f[istart])); end=int(float(f[iend]))
    except (ValueError,IndexError):
        continue
    name=f[iname] if 0<=iname<len(f) else f"sv_{chrom}_{start}"
    if isvlen>=0 and isvlen<len(f):
        try:
            if abs(int(float(f[isvlen]))) < '"$MIN_SV_LEN"': continue
        except ValueError: pass
    else:
        if end-start < '"$MIN_SV_LEN"': continue
    def af(i):
        if 0<=i<len(f):
            try: return float(f[i])
            except ValueError: return 0.0
        return 0.0
    diff=abs(af(iafr)-af(ieas))
    # normalise chrom to "chrN" (HPRC/cactus GRCh38 paths use chr-prefixed names)
    if not chrom.startswith("chr"): chrom="chr"+chrom
    out.write(f"{chrom}\t{start}\t{end}\t{name}\t{diff:.4f}\n")
' | sort -k1,1 -k2,2n > "$GNOMAD_BED"
    log "  wrote $(wc -l < "$GNOMAD_BED") gnomAD-SV records"
fi

# ── 3. ClinVar SVs (dbVar nstd102) -> BED (chrom start end id) ───────────────
log "=== 3/3: ClinVar SVs (dbVar nstd102) ==="
if [[ -s "$CLINVAR_BED" ]]; then
    log "  clinvar_sv.bed already present (skip)"
else
    fetch "$CLINVAR_URL" "$DL_DIR/nstd102.vcf.gz"
    log "  converting -> $CLINVAR_BED"
    # VCF -> BED. Use END from INFO when present (SVs), else POS+len(REF).
    zcat "$DL_DIR/nstd102.vcf.gz" | python3 -c '
import sys
ml='"$MIN_SV_LEN"'
for line in sys.stdin:
    if line.startswith("#"): continue
    f=line.rstrip("\n").split("\t")
    if len(f)<8: continue
    chrom,pos,vid,ref,alt=f[0],f[1],f[2],f[3],f[4]
    info=f[7]
    try: start=int(pos)-1  # VCF is 1-based; BED is 0-based
    except ValueError: continue
    end=None
    for kv in info.split(";"):
        if kv.startswith("END="):
            try: end=int(kv[4:])
            except ValueError: pass
    if end is None:
        end=start+max(len(ref),1)
    if end-start < ml: continue
    if not chrom.startswith("chr"): chrom="chr"+chrom
    name=vid if vid and vid!="." else f"clinvar_{chrom}_{start}"
    print(f"{chrom}\t{start}\t{end}\t{name}")
' | sort -k1,1 -k2,2n > "$CLINVAR_BED"
    log "  wrote $(wc -l < "$CLINVAR_BED") ClinVar SV records"
fi

log "=== reference data ready in $REF_DIR ==="
log "  GRCh38:    $GRCH38"
log "  gnomAD-SV: $GNOMAD_BED  (cols: chrom start end id fst=|AF_afr-AF_eas|)"
log "  ClinVar:   $CLINVAR_BED (cols: chrom start end id)"
log ""
log "exp6 reads these via G2G_GRCH38 / G2G_GNOMAD_SV / G2G_CLINVAR_SV"
log "(defaults already point at \$REF_DIR/{GRCh38.fa,gnomad_sv_v4.bed,clinvar_sv.bed})."
log ""
log "You can delete the raw downloads to save space: rm -rf $DL_DIR"
