#!/bin/bash --login
# =============================================================================
# Experiment 6 — Biological Validation   (Q5: Does G2G find real biology?)
# =============================================================================
# Roadmap §2.7.4 / §3.5 (Exp 6). This is the experiment that answers the
# question every systems reviewer asks of a pangenome paper: "why pangenomes
# and not some other irregular graph workload?" Without it, G2G is a fast
# graph engine looking for a problem. With it, G2G recovers known
# population-differentiated structural variation directly from a graph-to-graph
# alignment — something sequence-to-reference pipelines cannot do.
#
# THE PIPELINE (six stages, each resumable)
# -----------------------------------------
#   S1  Resolve + verify HPRC v2.0 assemblies for two population panels:
#         AFR  (n=15 samples -> 30 haplotypes)
#         EAS  (n=15 samples -> 30 haplotypes)
#       Assemblies are PRE-STAGED on the login node (prestage_hprc.sh) and
#       verified here by SHA-256 against the official HPRC index. Resumable.
#   S2  Build two pangenome graphs with minigraph-cactus (the HPRC-recommended
#         pipeline). ~48 CPU-hours per graph. Cached: skipped if GFA exists.
#   S3  Run G2G (parallel, 128 threads) to align GA(AFR) <-> GB(EAS).
#   S4  Lift delta alleles to GRCh38 coordinates (BED).
#   S5  Intersect with gnomAD-SV + ClinVar using bedtools.
#   S6  Permutation enrichment + a vg-giraffe delta-set comparison.
#
# WHY THIS IS A SEPARATE, LONG-RUNNING JOB
#   S2 (graph construction, ~96 CPU-h total) is the long pole. It is designed
#   to run unattended and resume after pre-emption.
#
# >>> NETWORK / DATA-STAGING NOTE (important — this is what was failing) <<<
#   The original Stage 1 hand-built S3 URLs and tried to wget them FROM THE
#   COMPUTE NODE. That failed for two independent reasons:
#     (1) the URLs were wrong — they omitted the
#         'year1_f1_assembly_v2_genbank/' directory, used the wrong base path
#         (an extra '/assemblies'), and used the non-genbank filename
#         'f1_assembly_v2.fa.gz' instead of 'f1_assembly_v2_genbank.fa.gz';
#     (2) ROCS testbed compute nodes (mcnode*) sit on a private VPN-gated
#         network with no verified outbound route to AWS S3, so each wget failed
#         INSTANTLY ("connection refused", ~1 s per file — exactly the symptom
#         in the exp6-*.err log).
#   FIX:
#     * Download all assemblies ONCE on the LOGIN NODE (mcmgt01) with
#         slurm/scripts/prestage_hprc.sh, which drives every URL from the
#         OFFICIAL HPRC index file (single source of truth) and SHA-verifies.
#     * This Stage 1 now only RESOLVES URLs/SHAs from that index and VERIFIES
#         the pre-staged files; it does not download from the compute node
#         unless G2G_ALLOW_NODE_DOWNLOAD=1 (only useful if this node truly has
#         S3 access).
#
# OUTPUTS  (all under results/exp6/)
#   delta_enrichment.tsv      per-truth-set fold-change + permutation p-value
#   overlap_counts.tsv        raw bedtools intersect counts (delta vs each set)
#   shared_pathogenic.tsv     ClinVar pathogenic SVs found in the SHARED graph
#   g2g_vs_giraffe.tsv        delta variants found by G2G but missed by giraffe
#   pangenome_stats.tsv       node/edge/snarl counts + build wall-time per panel
#   delta.bed / shared.bed    lifted intervals (reproducibility)
#   env.txt, manifest.sha256
#
# SUBMIT (after pre-staging on the login node):
#   bash   slurm/scripts/prestage_hprc.sh      # login node, one time
#   sbatch slurm/scripts/exp6_biological.sh    # compute node
# Long partition recommended (build stage can exceed 24h on its own).
# =============================================================================

#SBATCH --job-name=g2g_exp6_bio
#SBATCH --output=slurm/logs/exp6-%j.out
#SBATCH --error=slurm/logs/exp6-%j.err
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --constraint=cpu_amd_epyc_7763

set -euo pipefail

# Resolve repo root robustly. Under SLURM the batch script is copied to a
# spool dir (/var/lib/slurm/...), so BASH_SOURCE no longer points into the repo.
# Prefer $G2G_REPO (exported by submit_all.sh / the user); fall back to the
# script's own location only for in-place runs.
if [[ -n "${G2G_REPO:-}" ]]; then
    REPO="$G2G_REPO"
else
    _src="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$_src/../.." && pwd)"
fi
SCRIPT_DIR="$REPO/slurm/scripts"
G2G="${G2G_BIN:-$REPO/build/g2g}"
OUTDIR="$REPO/results/exp6"
NCPU="$(nproc)"

# ── Data roots (large; point at the >30 TB scratch, NOT $HOME) ───────────────
DATA_ROOT="${G2G_DATA_ROOT:-$REPO/data/hprc}"
ASM_DIR="$DATA_ROOT/assemblies"        # raw HPRC .fa.gz assemblies (pre-staged)
GRAPH_DIR="$DATA_ROOT/pangenomes"      # built GFA graphs
REF_DIR="$DATA_ROOT/reference"         # GRCh38 + truth sets
WORK_DIR="${G2G_SCRATCH:-$DATA_ROOT/work}"   # cactus scratch (huge, ephemeral)

# ── Truth-set locations (user supplies or pre-stages) ────────────────────────
GNOMAD_SV="${G2G_GNOMAD_SV:-$REF_DIR/gnomad_sv_v4.bed}"
CLINVAR_SV="${G2G_CLINVAR_SV:-$REF_DIR/clinvar_sv.bed}"
GRCH38="${G2G_GRCH38:-$REF_DIR/GRCh38.fa}"

# ── Panel definition ─────────────────────────────────────────────────────────
AFR_SAMPLES="${G2G_AFR_LIST:-$SCRIPT_DIR/panels/afr_samples.txt}"
EAS_SAMPLES="${G2G_EAS_LIST:-$SCRIPT_DIR/panels/eas_samples.txt}"

# ── Tunables ─────────────────────────────────────────────────────────────────
PERM_N="${G2G_PERM_N:-1000}"           # permutations for the null distribution
MIN_SV_LEN="${G2G_MIN_SV_LEN:-50}"     # SV size floor (bp); below this is not an SV
FST_HI="${G2G_FST_HI:-0.2}"            # high-F_ST threshold for gnomAD-SV stratum

mkdir -p "$OUTDIR" "$ASM_DIR" "$GRAPH_DIR" "$REF_DIR" "$WORK_DIR" "$REPO/slurm/logs"
source "$SCRIPT_DIR/lib_measure.sh"

# ── conda ────────────────────────────────────────────────────────────────────
if [[ -n "${CONDA_BASE:-}" && -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
fi
conda activate "${G2G_ENV:-phase5}" 2>/dev/null || log "WARN: conda env not activated"

capture_env "$OUTDIR/env.txt"
MANIFEST="$OUTDIR/manifest.sha256"; : > "$MANIFEST"
log "exp6 start  ncpu=$NCPU  data_root=$DATA_ROOT"

# Per-run result table: always start clean (header only), so re-running the job
# after a pre-emption does not accumulate duplicate rows.
fresh_tsv() { local path="$1"; shift; printf '%b\n' "$1" > "$path"; }

have() { command -v "$1" >/dev/null 2>&1; }
HAVE_CACTUS=$(have cactus-pangenome && echo 1 || echo 0)
HAVE_MINIGRAPH=$(have minigraph && echo 1 || echo 0)
HAVE_BEDTOOLS=$(have bedtools && echo 1 || echo 0)
HAVE_VG=$(have vg && echo 1 || echo 0)
HAVE_SAMTOOLS=$(have samtools && echo 1 || echo 0)
log "tools: cactus=$HAVE_CACTUS minigraph=$HAVE_MINIGRAPH bedtools=$HAVE_BEDTOOLS vg=$HAVE_VG samtools=$HAVE_SAMTOOLS"

# =============================================================================
# STAGE 1 — RESOLVE + VERIFY HPRC ASSEMBLIES  (index-driven; pre-stage aware)
# =============================================================================
# See the NETWORK / DATA-STAGING NOTE in the header. Files are pre-staged on the
# login node by prestage_hprc.sh; here we resolve their canonical URLs + SHA-256
# from the official HPRC index and verify the staged copies. We do NOT download
# from the compute node unless G2G_ALLOW_NODE_DOWNLOAD=1.
DLLOG="$OUTDIR/download_log.tsv"
fresh_tsv "$DLLOG" "panel\tsample\thaplotype\tfile\tbytes\tsha256\tstatus"

# Local copy of the HPRC index (pre-staged by prestage_hprc.sh). If absent and
# this node has internet (G2G_ALLOW_NODE_DOWNLOAD=1), fetch it; otherwise S1
# cannot resolve URLs and prints an actionable hint.
HPRC_INDEX_URL="${HPRC_INDEX_URL:-https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index}"
INDEX_FILE="${G2G_HPRC_INDEX:-$DATA_ROOT/index/Year1_assemblies_v2_genbank.index}"
ALLOW_NODE_DOWNLOAD="${G2G_ALLOW_NODE_DOWNLOAD:-0}"

mkdir -p "$(dirname "$INDEX_FILE")"
if [[ ! -s "$INDEX_FILE" && "$ALLOW_NODE_DOWNLOAD" -eq 1 ]]; then
    log "STAGE1: fetching HPRC index (node download allowed) -> $INDEX_FILE"
    wget --tries=5 --waitretry=10 --timeout=60 -q -O "$INDEX_FILE" "$HPRC_INDEX_URL" || true
fi

INDEX_OK=0
if [[ -s "$INDEX_FILE" ]] && head -1 "$INDEX_FILE" | grep -q 'hap1_aws_fasta'; then
    INDEX_OK=1
    log "STAGE1: HPRC index ready ($(wc -l < "$INDEX_FILE") rows) -> $INDEX_FILE"
else
    log "STAGE1: HPRC index NOT available at $INDEX_FILE"
    log "STAGE1:   Run 'bash slurm/scripts/prestage_hprc.sh' on the LOGIN NODE first,"
    log "STAGE1:   or set G2G_ALLOW_NODE_DOWNLOAD=1 if this node can reach AWS S3."
fi

# Convert an s3:// HPRC URL to the no-egress https endpoint.
s3_to_https() { sed 's%^s3://human-pangenomics/%https://s3-us-west-2.amazonaws.com/human-pangenomics/%'; }

# Resolve (url, sha) for one sample+haplotype from the index.
# paternal => hap1 (url col 2, sha col 6); maternal => hap2 (url col 3, sha col 7)
resolve_index() {
    local sample="$1" hap="$2" urlcol shacol
    case "$hap" in
        paternal) urlcol=2; shacol=6 ;;
        maternal) urlcol=3; shacol=7 ;;
        *) return 2 ;;
    esac
    awk -F'\t' -v s="$sample" -v uc="$urlcol" -v sc="$shacol" '
        $1==s { print $uc "\t" $sc; found=1 }
        END   { if(!found) exit 3 }' "$INDEX_FILE"
}

# Verify (and optionally download) one haplotype assembly. Echoes its path on success.
fetch_assembly() {
    local panel="$1" sample="$2" hap="$3"

    if [[ "$INDEX_OK" -eq 0 ]]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$panel" "$sample" "$hap" "NA" "0" "NA" "no_index" >> "$DLLOG"
        return 1
    fi

    local row s3url sha
    if ! row=$(resolve_index "$sample" "$hap"); then
        log "  [$panel] $sample.$hap NOT in HPRC index — check the sample ID"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$panel" "$sample" "$hap" "NA" "0" "NA" "not_in_index" >> "$DLLOG"
        return 1
    fi
    s3url=$(printf '%s' "$row" | cut -f1)
    sha=$(  printf '%s' "$row" | cut -f2)
    local url; url=$(printf '%s' "$s3url" | s3_to_https)
    local fname; fname=$(basename "$url")
    local dest="$ASM_DIR/$fname"

    # Already pre-staged + verified?
    if [[ -s "$dest" ]] && [[ "$(sha_of "$dest")" == "$sha" ]]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$panel" "$sample" "$hap" "$fname" "$(stat -c%s "$dest")" "$sha" "verified" >> "$DLLOG"
        log "  [$panel] verified  $fname"
        echo "$dest"; return 0
    fi

    # Present but failed verification.
    if [[ -s "$dest" ]]; then
        log "  [$panel] $fname present but SHA mismatch — re-stage it"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$panel" "$sample" "$hap" "$fname" "$(stat -c%s "$dest")" "$(sha_of "$dest")" "sha_mismatch" >> "$DLLOG"
        [[ "$ALLOW_NODE_DOWNLOAD" -eq 1 ]] || return 1
    fi

    # Missing. Only download from THIS node if explicitly allowed.
    if [[ "$ALLOW_NODE_DOWNLOAD" -ne 1 ]]; then
        log "  [$panel] MISSING $fname — pre-stage on the login node (prestage_hprc.sh)"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$panel" "$sample" "$hap" "$fname" "0" "$sha" "missing_prestage" >> "$DLLOG"
        return 1
    fi

    log "  [$panel] fetching $fname"
    if wget --continue --tries=5 --waitretry=30 --timeout=120 -q -O "$dest" "$url" 2>>"$OUTDIR/wget.err"; then
        local got bytes; got=$(sha_of "$dest"); bytes=$(stat -c%s "$dest" 2>/dev/null || echo 0)
        if [[ "$got" != "$sha" ]]; then
            log "  [$panel] SHA mismatch after download for $fname — deleting"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$panel" "$sample" "$hap" "$fname" "$bytes" "$got" "sha_mismatch" >> "$DLLOG"
            rm -f "$dest"; return 1
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$panel" "$sample" "$hap" "$fname" "$bytes" "$got" "downloaded" >> "$DLLOG"
        echo "$dest"; return 0
    else
        log "  [$panel] download FAILED for $fname (does this node have S3 access?)"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$panel" "$sample" "$hap" "$fname" "0" "NA" "failed" >> "$DLLOG"
        return 1
    fi
}

# Build a cactus seqfile (TSV: <name>\t<path>) for one panel from verified files.
build_seqfile() {
    local panel="$1" listfile="$2" seqfile="$3"
    : > "$seqfile"
    if [[ ! -s "$listfile" ]]; then
        log "STAGE1: panel list $listfile missing/empty — cannot build $panel"
        return 1
    fi
    # Reference row first (minigraph-cactus convention).
    if [[ -s "$GRCH38" ]]; then
        printf '_MINIGRAPH_\t%s\n' "$GRCH38" >> "$seqfile"
    fi
    local n_ok=0
    while IFS= read -r sample; do
        [[ -z "$sample" || "$sample" == \#* ]] && continue
        for hap in paternal maternal; do
            if dest=$(fetch_assembly "$panel" "$sample" "$hap"); then
                printf '%s.%s\t%s\n' "$sample" "$hap" "$dest" >> "$seqfile"
                echo "$(sha_of "$dest")  $(basename "$dest")" >> "$MANIFEST"
                n_ok=$((n_ok+1))
            fi
        done
    done < "$listfile"
    log "STAGE1: $panel — $n_ok haplotype assemblies ready -> $seqfile"
    [[ "$n_ok" -ge 2 ]]
}

log "=== STAGE 1: HPRC assembly resolve + verify ==="
AFR_SEQFILE="$DATA_ROOT/afr.seqfile"
EAS_SEQFILE="$DATA_ROOT/eas.seqfile"
S1_OK=1
build_seqfile "AFR" "$AFR_SAMPLES" "$AFR_SEQFILE" || S1_OK=0
build_seqfile "EAS" "$EAS_SAMPLES" "$EAS_SEQFILE" || S1_OK=0

if [[ "$S1_OK" -eq 0 ]]; then
    log "STAGE1: one or both panels incomplete."
    log "STAGE1: If assemblies are missing, pre-stage them on the LOGIN NODE:"
    log "STAGE1:     bash slurm/scripts/prestage_hprc.sh"
    log "STAGE1: then re-submit this job. (See download_log.tsv for per-file status.)"
fi

# =============================================================================
# STAGE 2 — BUILD PANGENOMES (minigraph-cactus; cached)
# =============================================================================
PGSTATS="$OUTDIR/pangenome_stats.tsv"
fresh_tsv "$PGSTATS" "panel\tgfa\tnodes\tedges\tpaths\tbuild_wall_sec\tpeak_rss_kb"

build_pangenome() {
    local panel="$1" seqfile="$2" out_gfa="$3"
    if [[ -s "$out_gfa" ]]; then
        log "STAGE2: $panel pangenome already built -> $out_gfa (skip)"
        return 0
    fi
    if [[ "$HAVE_CACTUS" -eq 0 && "$HAVE_MINIGRAPH" -eq 0 ]]; then
        log "STAGE2: neither cactus-pangenome nor minigraph available — cannot build $panel"
        return 1
    fi
    local cactus_out="$GRAPH_DIR/${panel,,}_cactus"
    mkdir -p "$cactus_out"
    log "STAGE2: building $panel pangenome (this is the ~48 CPU-h stage)"

    local m
    if [[ "$HAVE_CACTUS" -eq 1 ]]; then
        m=$(measure_full cactus-pangenome "$WORK_DIR/js_${panel,,}" "$seqfile" \
                --outDir "$cactus_out" --outName "${panel,,}" \
                --reference _MINIGRAPH_ --gfa --gbz --threads "$NCPU" \
                --mgCores "$NCPU" --mapCores "$NCPU")
        if [[ -s "$cactus_out/${panel,,}.gfa.gz" ]]; then
            gunzip -c "$cactus_out/${panel,,}.gfa.gz" > "$out_gfa"
        fi
    else
        local refs; refs=$(awk '{print $2}' "$seqfile")
        m=$(measure_full bash -c "minigraph -cxggs -t $NCPU $refs > '$out_gfa'")
    fi
    read -r wall peak user sys cpu ec <<< "$m"

    if [[ ! -s "$out_gfa" ]]; then
        log "STAGE2: $panel build produced no GFA (exit=$ec)"; return 1
    fi
    local nodes edges paths
    nodes=$(grep -c '^S' "$out_gfa" 2>/dev/null || echo 0)
    edges=$(grep -c '^L' "$out_gfa" 2>/dev/null || echo 0)
    paths=$(grep -cE '^(P|W)' "$out_gfa" 2>/dev/null || echo 0)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$panel" "$out_gfa" "$nodes" "$edges" "$paths" "$wall" "$peak" >> "$PGSTATS"
    echo "$(sha_of "$out_gfa")  $(basename "$out_gfa")" >> "$MANIFEST"
    log "STAGE2: $panel done  nodes=$nodes edges=$edges paths=$paths wall=${wall}s"
}

log "=== STAGE 2: pangenome construction ==="
AFR_GFA="$GRAPH_DIR/afr.gfa"
EAS_GFA="$GRAPH_DIR/eas.gfa"
S2_OK=1
if [[ "$S1_OK" -eq 1 ]]; then
    build_pangenome "AFR" "$AFR_SEQFILE" "$AFR_GFA" || S2_OK=0
    build_pangenome "EAS" "$EAS_SEQFILE" "$EAS_GFA" || S2_OK=0
else
    log "STAGE2: skipped — Stage 1 did not complete both panels"
    S2_OK=0
fi

# =============================================================================
# STAGE 3 — RUN G2G (parallel) AFR <-> EAS
# =============================================================================
G2GLOG="$OUTDIR/g2g_run.tsv"
fresh_tsv "$G2GLOG" "stage\twall_sec\tpeak_rss_kb\tcpu_pct\texit"
DELTA_GFA="$OUTDIR/delta.gfa"
SHARED_GFA="$OUTDIR/shared.gfa"
ALIGN_GAF="$OUTDIR/afr_eas.gaf"

log "=== STAGE 3: G2G alignment AFR<->EAS ==="
S3_OK=1
if [[ "$S2_OK" -eq 1 && -x "$G2G" ]]; then
    stderr_f="$OUTDIR/g2g_run.stderr"
    m=$(measure_full_capture "$stderr_f" "$G2G" "$AFR_GFA" "$EAS_GFA" \
            --gaf "$ALIGN_GAF" --delta-gfa "$DELTA_GFA" --shared-gfa "$SHARED_GFA" \
            --threads "$NCPU")
    read -r wall peak user sys cpu ec <<< "$m"
    printf 'align\t%s\t%s\t%s\t%s\n' "$wall" "$peak" "$cpu" "$ec" >> "$G2GLOG"
    log "STAGE3: g2g exit=$ec wall=${wall}s rss=${peak}kB"
    if grep -q '\[timing\]' "$stderr_f" 2>/dev/null; then
        grep '\[timing\]\|\[runstats\]' "$stderr_f" | while read -r ln; do log "  $ln"; done
    fi
    [[ "$ec" == "0" && -s "$DELTA_GFA" ]] || S3_OK=0
    for g in "$DELTA_GFA" "$SHARED_GFA"; do
        [[ -s "$g" ]] && echo "$(sha_of "$g")  $(basename "$g")" >> "$MANIFEST"
    done
else
    log "STAGE3: skipped — pangenomes unavailable or G2G binary missing at $G2G"
    S3_OK=0
fi

# =============================================================================
# STAGE 4 — LIFT DELTA/SHARED ALLELES TO GRCh38 COORDINATES (BED)
# =============================================================================
DELTA_BED="$OUTDIR/delta.bed"
SHARED_BED="$OUTDIR/shared.bed"

lift_to_bed() {
    local in_gfa="$1" out_bed="$2" label="$3"
    : > "$out_bed"
    [[ -s "$in_gfa" ]] || { log "STAGE4: $label graph empty — no BED"; return 1; }

    if "$G2G" lift --help >/dev/null 2>&1; then
        "$G2G" lift --gfa "$in_gfa" --ref-prefix GRCh38 \
               --min-len "$MIN_SV_LEN" --out "$out_bed" 2>>"$OUTDIR/lift.err" || true
    fi
    if [[ ! -s "$out_bed" ]]; then
        awk -v ml="$MIN_SV_LEN" '
            $1=="W" && $2 ~ /GRCh38/ {
                start=$5+0; end=$6+0;
                if (end-start >= ml) printf "%s\t%d\t%d\t%s\n", $4, start, end, "delta";
            }' "$in_gfa" | sort -k1,1 -k2,2n > "$out_bed" 2>/dev/null || true
    fi
    local n; n=$(wc -l < "$out_bed" 2>/dev/null || echo 0)
    log "STAGE4: $label -> $out_bed ($n intervals >= ${MIN_SV_LEN}bp)"
    echo "$(sha_of "$out_bed")  $(basename "$out_bed")" >> "$MANIFEST"
    [[ "$n" -gt 0 ]]
}

log "=== STAGE 4: lift alleles to GRCh38 BED ==="
S4_OK=1
if [[ "$S3_OK" -eq 1 ]]; then
    lift_to_bed "$DELTA_GFA"  "$DELTA_BED"  "delta"  || S4_OK=0
    lift_to_bed "$SHARED_GFA" "$SHARED_BED" "shared" || true
else
    log "STAGE4: skipped — no delta graph"; S4_OK=0
fi

# =============================================================================
# STAGE 5 — INTERSECT WITH TRUTH SETS (bedtools)
# =============================================================================
OVERLAPS="$OUTDIR/overlap_counts.tsv"
fresh_tsv "$OVERLAPS" "query_set\ttruth_set\tquery_n\toverlap_n\tfraction"
SHARED_PATH="$OUTDIR/shared_pathogenic.tsv"

intersect_count() {
    local query_bed="$1" truth_bed="$2"
    local qn on
    qn=$(wc -l < "$query_bed" 2>/dev/null || echo 0)
    if [[ "$HAVE_BEDTOOLS" -eq 1 && -s "$truth_bed" && "$qn" -gt 0 ]]; then
        on=$(bedtools intersect -u -a "$query_bed" -b "$truth_bed" 2>/dev/null | wc -l)
    else
        on=0
    fi
    echo "$qn $on"
}

record_overlap() {
    local qset="$1" tset="$2" qbed="$3" tbed="$4"
    read -r qn on <<< "$(intersect_count "$qbed" "$tbed")"
    local frac="0"
    [[ "$qn" -gt 0 ]] && frac=$(awk "BEGIN{printf \"%.4f\", $on/$qn}")
    printf '%s\t%s\t%s\t%s\t%s\n' "$qset" "$tset" "$qn" "$on" "$frac" >> "$OVERLAPS"
    log "STAGE5: $qset ∩ $tset -> $on/$qn ($frac)"
}

log "=== STAGE 5: truth-set intersection ==="
S5_OK=1
if [[ "$S4_OK" -eq 1 ]]; then
    GNOMAD_HI="$OUTDIR/gnomad_sv_highfst.bed"
    if [[ -s "$GNOMAD_SV" ]]; then
        awk -v t="$FST_HI" 'NF>=5 && ($5+0)>=t {print $1"\t"$2"\t"$3"\t"$4}' \
            "$GNOMAD_SV" > "$GNOMAD_HI" 2>/dev/null || : > "$GNOMAD_HI"
    fi

    record_overlap "delta"  "gnomad_sv_all"    "$DELTA_BED" "$GNOMAD_SV"
    [[ -s "$GNOMAD_HI" ]] && record_overlap "delta" "gnomad_sv_highfst" "$DELTA_BED" "$GNOMAD_HI"
    record_overlap "delta"  "clinvar"          "$DELTA_BED" "$CLINVAR_SV"
    [[ -s "$SHARED_BED" ]] && record_overlap "shared" "clinvar" "$SHARED_BED" "$CLINVAR_SV"

    if [[ "$HAVE_BEDTOOLS" -eq 1 && -s "$SHARED_BED" && -s "$CLINVAR_SV" ]]; then
        fresh_tsv "$SHARED_PATH" "chrom\tstart\tend\tclinvar_id"
        bedtools intersect -u -a "$CLINVAR_SV" -b "$SHARED_BED" 2>/dev/null \
            | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' >> "$SHARED_PATH" || true
        log "STAGE5: shared∩ClinVar pathogenic dumped -> $SHARED_PATH"
    fi
else
    log "STAGE5: skipped — no lifted intervals"; S5_OK=0
fi

# =============================================================================
# STAGE 6a — PERMUTATION ENRICHMENT TEST
# =============================================================================
ENRICH="$OUTDIR/delta_enrichment.tsv"
fresh_tsv "$ENRICH" "truth_set\tobserved\tnull_mean\tnull_sd\tfold_change\tp_value\tperms"

permutation_test() {
    local query_bed="$1" truth_bed="$2" tset="$3"
    if [[ "$HAVE_BEDTOOLS" -eq 0 || ! -s "$query_bed" || ! -s "$truth_bed" ]]; then
        log "STAGE6: permutation skipped for $tset (inputs/tools missing)"
        return 0
    fi
    local genome="$REF_DIR/GRCh38.genome"
    if [[ ! -s "$genome" ]]; then
        if [[ "$HAVE_SAMTOOLS" -eq 1 && -s "$GRCH38" ]]; then
            [[ -s "${GRCH38}.fai" ]] || samtools faidx "$GRCH38" 2>/dev/null || true
            cut -f1,2 "${GRCH38}.fai" > "$genome" 2>/dev/null || true
        fi
    fi
    if [[ ! -s "$genome" ]]; then
        log "STAGE6: no genome sizes file — cannot shuffle; reporting observed only"
        local obs; read -r _ obs <<< "$(intersect_count "$query_bed" "$truth_bed")"
        printf '%s\t%s\tNA\tNA\tNA\tNA\t0\n' "$tset" "$obs" >> "$ENRICH"
        return 0
    fi

    local observed; read -r _ observed <<< "$(intersect_count "$query_bed" "$truth_bed")"
    log "STAGE6: $tset observed overlap = $observed; running $PERM_N permutations"

    local nulls="$OUTDIR/.null_${tset}.txt"; : > "$nulls"
    local incl=()
    [[ -s "$REF_DIR/mappable.bed" ]] && incl=(-incl "$REF_DIR/mappable.bed")
    for ((i=0; i<PERM_N; i++)); do
        bedtools shuffle -i "$query_bed" -g "$genome" -noOverlapping "${incl[@]}" 2>/dev/null \
            | bedtools intersect -u -a - -b "$truth_bed" 2>/dev/null | wc -l >> "$nulls"
    done

    python3 - "$nulls" "$observed" "$tset" "$PERM_N" "$ENRICH" <<'PYEOF'
import sys, statistics as st
nulls_f, observed, tset, perms, out = sys.argv[1:6]
observed = float(observed); perms = int(perms)
vals = [float(x) for x in open(nulls_f) if x.strip()]
if not vals:
    sys.exit(0)
mean = st.mean(vals)
sd = st.pstdev(vals) if len(vals) > 1 else 0.0
ge = sum(1 for v in vals if v >= observed)
p = (ge + 1) / (len(vals) + 1)
fold = observed / mean if mean > 0 else float("inf")
with open(out, "a") as f:
    f.write(f"{tset}\t{observed:.0f}\t{mean:.2f}\t{sd:.2f}\t{fold:.3f}\t{p:.5f}\t{len(vals)}\n")
print(f"[perm] {tset}: obs={observed:.0f} null={mean:.2f}±{sd:.2f} fold={fold:.2f} p={p:.4f}")
PYEOF
    rm -f "$nulls"
}

log "=== STAGE 6a: permutation enrichment ==="
if [[ "$S5_OK" -eq 1 ]]; then
    GNOMAD_HI="$OUTDIR/gnomad_sv_highfst.bed"
    [[ -s "$GNOMAD_HI" ]] && permutation_test "$DELTA_BED" "$GNOMAD_HI" "gnomad_sv_highfst"
    permutation_test "$DELTA_BED" "$GNOMAD_SV"  "gnomad_sv_all"
    permutation_test "$DELTA_BED" "$CLINVAR_SV" "clinvar"
else
    log "STAGE6a: skipped"
fi

# =============================================================================
# STAGE 6b — G2G vs vg-giraffe DELTA-SET COMPARISON
# =============================================================================
GVG="$OUTDIR/g2g_vs_giraffe.tsv"
fresh_tsv "$GVG" "metric\tvalue"

log "=== STAGE 6b: G2G vs vg-giraffe ==="
if [[ "$S4_OK" -eq 1 && "$HAVE_VG" -eq 1 && "$HAVE_BEDTOOLS" -eq 1 && -s "$AFR_GFA" ]]; then
    GIRAFFE_BED="$OUTDIR/giraffe_found.bed"
    stderr_f="$OUTDIR/giraffe.stderr"
    m=$(measure_full_capture "$stderr_f" bash -c "
        set -e
        vg autoindex --workflow giraffe -g '$AFR_GFA' -p '$GRAPH_DIR/afr_idx' -t $NCPU 2>/dev/null
        : > '$OUTDIR/eas_haps.fa'
        awk '\$1==\"W\" && \$2 !~ /GRCh38/ {next}' /dev/null 2>/dev/null || true
        vg paths -x '$EAS_GFA' -F 2>/dev/null > '$OUTDIR/eas_haps.fa' || true
        vg giraffe -Z '$GRAPH_DIR/afr_idx.giraffe.gbz' -f '$OUTDIR/eas_haps.fa' -t $NCPU 2>/dev/null \
            | vg surject -x '$GRAPH_DIR/afr_idx.giraffe.gbz' -b - 2>/dev/null \
            | bedtools bamtobed -i - 2>/dev/null \
            | sort -k1,1 -k2,2n > '$GIRAFFE_BED' || true
    ")
    read -r wall peak user sys cpu ec <<< "$m"

    if [[ -s "$GIRAFFE_BED" && -s "$DELTA_BED" ]]; then
        g2g_only=$(bedtools intersect -v -a "$DELTA_BED" -b "$GIRAFFE_BED" 2>/dev/null | wc -l)
        delta_n=$(wc -l < "$DELTA_BED")
        giraffe_n=$(wc -l < "$GIRAFFE_BED")
        frac="0"; [[ "$delta_n" -gt 0 ]] && frac=$(awk "BEGIN{printf \"%.4f\", $g2g_only/$delta_n}")
        {
            printf 'delta_intervals\t%s\n' "$delta_n"
            printf 'giraffe_intervals\t%s\n' "$giraffe_n"
            printf 'g2g_unique\t%s\n' "$g2g_only"
            printf 'g2g_unique_fraction\t%s\n' "$frac"
            printf 'giraffe_wall_sec\t%s\n' "$wall"
        } >> "$GVG"
        log "STAGE6b: G2G-unique delta intervals = $g2g_only/$delta_n ($frac); giraffe wall=${wall}s"
    else
        log "STAGE6b: giraffe produced no intervals — comparison inconclusive"
        printf 'status\tinconclusive\n' >> "$GVG"
    fi
else
    log "STAGE6b: skipped (vg/bedtools/graph unavailable)"
    printf 'status\tskipped\n' >> "$GVG"
fi

# =============================================================================
# DONE — summary
# =============================================================================
log "=== Exp6 DONE ==="
log "  download log:     $DLLOG"
log "  pangenome stats:  $PGSTATS"
log "  g2g run:          $G2GLOG"
log "  delta BED:        $DELTA_BED"
log "  overlaps:         $OVERLAPS"
log "  enrichment:       $ENRICH"
log "  g2g vs giraffe:   $GVG"
log "  shared pathogenic:$SHARED_PATH"
log "  manifest:         $MANIFEST"
log "stage status: S1=$S1_OK S2=$S2_OK S3=$S3_OK S4=$S4_OK S5=$S5_OK"

if [[ ! -s "$ENRICH" || $(wc -l < "$ENRICH") -le 1 ]]; then
    log "WARN: no enrichment rows produced — biological validation incomplete"
    exit 3
fi
log "exp6 completed with enrichment results."
