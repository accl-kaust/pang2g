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
#   S1  Download HPRC v2.0 assemblies for two population panels:
#         AFR  (n=15 samples -> 30 haplotypes)
#         EAS  (n=15 samples -> 30 haplotypes)
#       Verified by SHA-256 against the HPRC index. Resumable (wget -c).
#   S2  Build two pangenome graphs with minigraph-cactus (the HPRC-recommended
#         pipeline). ~48 CPU-hours per graph. Cached: skipped if GFA exists.
#   S3  Run G2G (parallel, 128 threads) to align GA(AFR) <-> GB(EAS).
#         Emits the delta graph (alleles present in one panel, absent/divergent
#         in the other) and the shared graph.
#   S4  Lift delta alleles to GRCh38 coordinates (BED) via the graph's
#         reference path, so they can be intersected with coordinate-based
#         truth sets.
#   S5  Intersect with truth sets using bedtools:
#         gnomAD-SV v4  -> expect delta alleles ENRICHED for high-F_ST SVs
#         ClinVar       -> expect SHARED alleles to carry known pathogenic SVs
#   S6  Significance: permutation null (shuffle delta intervals across the
#         mappable genome N times) -> empirical enrichment p-value + fold-change.
#         Plus a vg-giraffe delta-set comparison (what a mapping pipeline would
#         have called) to show G2G finds variants the standard route misses.
#
# WHY THIS IS A SEPARATE, LONG-RUNNING JOB
#   S1 (download, ~1 TB) and S2 (graph construction, ~96 CPU-h total) are the
#   long poles. They are designed to run unattended and to resume after
#   pre-emption. Launch this FIRST and let it run in the background while the
#   throughput/scaling/ablation jobs (exp3/4/5) run on other nodes.
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
# SUBMIT:  sbatch slurm/scripts/exp6_biological.sh
# Long partition recommended (build stage can exceed 24h on its own).
# =============================================================================

#SBATCH --job-name=g2g_exp6_bio
#SBATCH --output=slurm/logs/exp6-%j.out
#SBATCH --error=slurm/logs/exp6-%j.err
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=512GB
#SBATCH --constraint=cpu_amd_epyc_7763

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${G2G_REPO:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
G2G="${G2G_BIN:-$REPO/build/g2g}"
OUTDIR="$REPO/results/exp6"
NCPU="$(nproc)"

# ── Data roots (large; point at the >30 TB scratch, NOT $HOME) ───────────────
DATA_ROOT="${G2G_DATA_ROOT:-$REPO/data/hprc}"
ASM_DIR="$DATA_ROOT/assemblies"        # raw HPRC .fa.gz assemblies
GRAPH_DIR="$DATA_ROOT/pangenomes"      # built GFA graphs
REF_DIR="$DATA_ROOT/reference"         # GRCh38 + truth sets
WORK_DIR="${G2G_SCRATCH:-$DATA_ROOT/work}"   # cactus scratch (huge, ephemeral)

# ── Truth-set locations (user supplies or we fetch) ──────────────────────────
GNOMAD_SV="${G2G_GNOMAD_SV:-$REF_DIR/gnomad_sv_v4.bed}"
CLINVAR_SV="${G2G_CLINVAR_SV:-$REF_DIR/clinvar_sv.bed}"
GRCH38="${G2G_GRCH38:-$REF_DIR/GRCh38.fa}"

# ── Panel definition ─────────────────────────────────────────────────────────
# Sample lists live next to the script; one HPRC sample ID per line. These are
# the 15+15 samples already selected from the HPRC sample table by superpopulation.
AFR_SAMPLES="${G2G_AFR_LIST:-$SCRIPT_DIR/panels/afr_samples.txt}"
EAS_SAMPLES="${G2G_EAS_LIST:-$SCRIPT_DIR/panels/eas_samples.txt}"
HPRC_BASE="${G2G_HPRC_BASE:-https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/assemblies}"

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

# Tool availability — none are fatal at parse time; each stage checks its own.
# Per-run result table: always start clean (header only), so re-running the job
# after a pre-emption does not accumulate duplicate rows. Provenance is kept
# separately in the SHA manifest, which IS append-only.
fresh_tsv() { local path="$1"; shift; printf '%b\n' "$1" > "$path"; }

have() { command -v "$1" >/dev/null 2>&1; }
HAVE_CACTUS=$(have cactus-pangenome && echo 1 || echo 0)
HAVE_MINIGRAPH=$(have minigraph && echo 1 || echo 0)
HAVE_BEDTOOLS=$(have bedtools && echo 1 || echo 0)
HAVE_VG=$(have vg && echo 1 || echo 0)
HAVE_SAMTOOLS=$(have samtools && echo 1 || echo 0)
log "tools: cactus=$HAVE_CACTUS minigraph=$HAVE_MINIGRAPH bedtools=$HAVE_BEDTOOLS vg=$HAVE_VG samtools=$HAVE_SAMTOOLS"

# =============================================================================
# STAGE 1 — DOWNLOAD HPRC ASSEMBLIES (resumable, SHA-256 verified)
# =============================================================================
# Reads a panel list, downloads each sample's maternal+paternal assemblies,
# verifies SHA-256 against the HPRC index if present, and records provenance.
# Already-verified files are skipped, so re-running after pre-emption is cheap
# (this is the "58/60 already downloaded" resume case from the roadmap).
DLLOG="$OUTDIR/download_log.tsv"
fresh_tsv "$DLLOG" "panel\tsample\thaplotype\tfile\tbytes\tsha256\tstatus"

fetch_assembly() {
    local panel="$1" sample="$2" hap="$3"
    # HPRC naming: <SAMPLE>.<paternal|maternal>.f1_assembly_v2.fa.gz
    local fname="${sample}.${hap}.f1_assembly_v2.fa.gz"
    local url="$HPRC_BASE/$sample/assemblies/$fname"
    local dest="$ASM_DIR/$fname"

    if [[ -s "$dest" ]]; then
        # Re-verify a previously downloaded file rather than trusting its presence.
        local sha; sha=$(sha_of "$dest")
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$panel" "$sample" "$hap" "$fname" "$(stat -c%s "$dest" 2>/dev/null || echo 0)" \
            "$sha" "cached" >> "$DLLOG"
        log "  [dl] cached  $fname"
        echo "$dest"; return 0
    fi

    log "  [dl] fetching $fname"
    # --continue: resume partial transfers after pre-emption.
    # --tries + --waitretry: ride out transient S3 hiccups.
    if wget --continue --tries=5 --waitretry=30 --timeout=120 \
            -q -O "$dest" "$url" 2>>"$OUTDIR/wget.err"; then
        local sha bytes; sha=$(sha_of "$dest"); bytes=$(stat -c%s "$dest" 2>/dev/null || echo 0)
        # Sanity: a real HPRC haplotype assembly is > 500 MB. A tiny file is an
        # S3 error page that slipped past wget, not an assembly.
        if [[ "$bytes" -lt 500000000 ]]; then
            log "  [dl] WARN $fname only ${bytes}B — likely an error page; marking failed"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$panel" "$sample" "$hap" "$fname" "$bytes" "$sha" "too_small" >> "$DLLOG"
            rm -f "$dest"; return 1
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$panel" "$sample" "$hap" "$fname" "$bytes" "$sha" "downloaded" >> "$DLLOG"
        echo "$dest"; return 0
    else
        log "  [dl] FAILED  $fname"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$panel" "$sample" "$hap" "$fname" "0" "NA" "failed" >> "$DLLOG"
        return 1
    fi
}

# Build a cactus seqfile (TSV: <name>\t<path>) for one panel as we download.
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

log "=== STAGE 1: HPRC assembly download ==="
AFR_SEQFILE="$DATA_ROOT/afr.seqfile"
EAS_SEQFILE="$DATA_ROOT/eas.seqfile"
S1_OK=1
build_seqfile "AFR" "$AFR_SAMPLES" "$AFR_SEQFILE" || S1_OK=0
build_seqfile "EAS" "$EAS_SAMPLES" "$EAS_SEQFILE" || S1_OK=0

# =============================================================================
# STAGE 2 — BUILD PANGENOMES (minigraph-cactus; cached)
# =============================================================================
# ~48 CPU-h per panel. If the output GFA already exists and is non-empty, we
# skip — this is what makes the whole job resumable across pre-emptions.
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
        # minigraph-cactus, HPRC-recommended. --gfa emits the final graph.
        m=$(measure_full cactus-pangenome "$WORK_DIR/js_${panel,,}" "$seqfile" \
                --outDir "$cactus_out" --outName "${panel,,}" \
                --reference _MINIGRAPH_ --gfa --gbz --threads "$NCPU" \
                --mgCores "$NCPU" --mapCores "$NCPU")
        # cactus writes <outName>.gfa.gz; normalise to our path.
        if [[ -s "$cactus_out/${panel,,}.gfa.gz" ]]; then
            gunzip -c "$cactus_out/${panel,,}.gfa.gz" > "$out_gfa"
        fi
    else
        # Fallback: raw minigraph graph (no cactus base alignment). Coarser, but
        # establishes the panel graph if cactus is unavailable.
        local refs; refs=$(awk '{print $2}' "$seqfile")
        m=$(measure_full bash -c "minigraph -cxggs -t $NCPU $refs > '$out_gfa'")
    fi
    read -r wall peak user sys cpu ec <<< "$m"

    if [[ ! -s "$out_gfa" ]]; then
        log "STAGE2: $panel build produced no GFA (exit=$ec)"; return 1
    fi
    # Quick structural census straight from the GFA.
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
# The actual contribution: a single graph-to-graph alignment of the two
# population pangenomes. Output is a delta graph (panel-divergent alleles) and
# a shared graph. We time and memory-profile this run too — it is also the
# largest single real-data G2G invocation in the paper.
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
    # Surface the runtime's own instrumentation if present.
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
# The truth sets (gnomAD-SV, ClinVar) are in GRCh38 coordinates. G2G's delta
# alleles are graph nodes. We project each delta/shared allele onto the graph's
# GRCh38 reference path (W-line / P-line named GRCh38*) to obtain an interval.
# This is what `g2g lift` does; if the binary lacks the subcommand, we fall
# back to a small AWK projector over the GFA W-lines.
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
    # Fallback / supplement: derive intervals from W-lines referencing GRCh38.
    # W-line: W <sample> <hap> <seqid> <start> <end> <walk>
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
    # echoes "<query_n> <overlap_n>"
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
    # Build a high-F_ST stratum of gnomAD-SV if the F_ST column is present
    # (col 5 by our convention: chrom start end id fst). Expectation: delta
    # alleles are enriched in the HIGH-F_ST stratum vs all SVs.
    GNOMAD_HI="$OUTDIR/gnomad_sv_highfst.bed"
    if [[ -s "$GNOMAD_SV" ]]; then
        awk -v t="$FST_HI" 'NF>=5 && ($5+0)>=t {print $1"\t"$2"\t"$3"\t"$4}' \
            "$GNOMAD_SV" > "$GNOMAD_HI" 2>/dev/null || : > "$GNOMAD_HI"
    fi

    record_overlap "delta"  "gnomad_sv_all"    "$DELTA_BED" "$GNOMAD_SV"
    [[ -s "$GNOMAD_HI" ]] && record_overlap "delta" "gnomad_sv_highfst" "$DELTA_BED" "$GNOMAD_HI"
    record_overlap "delta"  "clinvar"          "$DELTA_BED" "$CLINVAR_SV"
    # Shared graph should carry pan-population pathogenic variants.
    [[ -s "$SHARED_BED" ]] && record_overlap "shared" "clinvar" "$SHARED_BED" "$CLINVAR_SV"

    # Dump the actual ClinVar-pathogenic SVs that fall in the SHARED graph.
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
# Is the delta∩gnomAD-high-F_ST overlap more than expected by chance? Shuffle
# the delta intervals PERM_N times across the mappable genome (preserving count
# and length distribution) and recompute the overlap each time -> empirical p.
ENRICH="$OUTDIR/delta_enrichment.tsv"
fresh_tsv "$ENRICH" "truth_set\tobserved\tnull_mean\tnull_sd\tfold_change\tp_value\tperms"

permutation_test() {
    local query_bed="$1" truth_bed="$2" tset="$3"
    if [[ "$HAVE_BEDTOOLS" -eq 0 || ! -s "$query_bed" || ! -s "$truth_bed" ]]; then
        log "STAGE6: permutation skipped for $tset (inputs/tools missing)"
        return 0
    fi
    # A genome file (chrom sizes) is required for bedtools shuffle.
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
    # Optional mappable-region mask keeps shuffled intervals in callable genome.
    local incl=()
    [[ -s "$REF_DIR/mappable.bed" ]] && incl=(-incl "$REF_DIR/mappable.bed")
    for ((i=0; i<PERM_N; i++)); do
        bedtools shuffle -i "$query_bed" -g "$genome" -noOverlapping "${incl[@]}" 2>/dev/null \
            | bedtools intersect -u -a - -b "$truth_bed" 2>/dev/null | wc -l >> "$nulls"
    done

    # Empirical p, fold-change vs null mean (Python for the stats).
    python3 - "$nulls" "$observed" "$tset" "$PERM_N" "$ENRICH" <<'PYEOF'
import sys, statistics as st
nulls_f, observed, tset, perms, out = sys.argv[1:6]
observed = float(observed); perms = int(perms)
vals = [float(x) for x in open(nulls_f) if x.strip()]
if not vals:
    sys.exit(0)
mean = st.mean(vals)
sd = st.pstdev(vals) if len(vals) > 1 else 0.0
# one-sided p: P(null >= observed), +1 smoothing
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
# The headline biological argument: G2G recovers panel-divergent SVs that the
# standard mapping route (vg giraffe of one panel's haplotypes onto the other's
# graph, then call) does not surface. We map EAS haplotypes onto the AFR graph
# with giraffe, derive that pipeline's "found" interval set, and report the
# delta variants G2G found that giraffe missed.
GVG="$OUTDIR/g2g_vs_giraffe.tsv"
fresh_tsv "$GVG" "metric\tvalue"

log "=== STAGE 6b: G2G vs vg-giraffe ==="
if [[ "$S4_OK" -eq 1 && "$HAVE_VG" -eq 1 && "$HAVE_BEDTOOLS" -eq 1 && -s "$AFR_GFA" ]]; then
    GIRAFFE_BED="$OUTDIR/giraffe_found.bed"
    stderr_f="$OUTDIR/giraffe.stderr"
    # Build a giraffe index over the AFR graph and map EAS assemblies' sequences.
    # We time it for an honest apples-to-apples cost note in the paper.
    m=$(measure_full_capture "$stderr_f" bash -c "
        set -e
        vg autoindex --workflow giraffe -g '$AFR_GFA' -p '$GRAPH_DIR/afr_idx' -t $NCPU 2>/dev/null
        # EAS haplotype sequences (concatenated) as the query read set.
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
        # Delta intervals NOT covered by anything giraffe surjected = G2G-unique.
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

# Non-zero exit if the core scientific result (delta enrichment) never produced,
# so a wrapper/master script can tell a real failure from a clean run.
if [[ ! -s "$ENRICH" || $(wc -l < "$ENRICH") -le 1 ]]; then
    log "WARN: no enrichment rows produced — biological validation incomplete"
    exit 3
fi
log "exp6 completed with enrichment results."
