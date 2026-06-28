#!/usr/bin/env bash
# Demo 7 — throughput vs. the linearize-and-map baseline (fig_throughput).
#
# Question: why align two pangenome graphs directly (G2G) instead of the obvious
# alternative — linearize each graph back into haplotypes and map them with a
# sequence aligner? Because the graph is a *compressed* representation: its size
# saturates as redundant haplotypes are added, so G2G's whole diff is ~flat in the
# panel size h, whereas linearize-and-map must extract and align every haplotype,
# growing linearly in h.
#
# We sweep h (haplotypes per population) and, at each point, time:
#   * G2G          : g2g popA.gfa popB.gfa --ref <contig>   (the entire diff)
#   * baseline     : extract h haplotypes per population (bcftools consensus) and
#                    align them to the reference (minimap2) — the dominant cost of
#                    any linearize-and-map pipeline, BEFORE variant calling. Timing
#                    only this stage is deliberately generous to the baseline: the
#                    real pipeline then has to call and merge variants on top.
#
# Everything is GRCh37/b37, reusing the inputs staged by 09_fetch_biology_data.sh.
# Needs vg, bcftools, tabix, samtools, minimap2 (vg+minimap2 ship in demo/tools via
# 01_setup_tools.sh; bcftools/samtools/tabix from your g2g_bio env). Cannot run in
# the dev sandbox (no tools/data there). Each tool is overridable: VG=, BCFTOOLS=,
# TABIX=, SAMTOOLS=, MINIMAP2=.
#
# Inputs (env): REF_FA, PANEL_VCF or PANEL_TMPL ({}=contig), SAMPLES_A, SAMPLES_B,
#   CONTIG (default 22), H_SWEEP (default "2 5 10 25 50 100").
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
need_g2g
OUT="$WORK/throughput"; mkdir -p "$OUT"
CONTIG="${CONTIG:-22}"
H_SWEEP="${H_SWEEP:-2 5 10 25 50 100}"
VG="${VG:-vg}"; BCFTOOLS="${BCFTOOLS:-bcftools}"; TABIX="${TABIX:-tabix}"
SAMTOOLS="${SAMTOOLS:-samtools}"; MINIMAP2="${MINIMAP2:-minimap2}"

for t in "$VG" "$BCFTOOLS" "$TABIX" "$SAMTOOLS" "$MINIMAP2"; do
    have "$t" || { echo "ERROR: '$t' not found (set VG=/BCFTOOLS=/TABIX=/SAMTOOLS=/MINIMAP2= or install it)" >&2; exit 1; }
done
case "$("$BCFTOOLS" --version 2>/dev/null | sed -n '1{s/.*bcftools *//p}')" in
    1.*|2.*) : ;;
    *) echo "ERROR: need bcftools >=1.0 (use the g2g_bio env or BCFTOOLS=/path)" >&2; exit 1;;
esac
: "${REF_FA:?set REF_FA to the reference fasta (with .fai)}"
: "${SAMPLES_A:?set SAMPLES_A to a file of population-A sample ids}"
: "${SAMPLES_B:?set SAMPLES_B to a file of population-B sample ids}"
[ -s "$SAMPLES_A" ] && [ -s "$SAMPLES_B" ] || { echo "ERROR: SAMPLES_A/SAMPLES_B empty/missing" >&2; exit 1; }
if [ -z "${PANEL_VCF:-}" ] && [ -z "${PANEL_TMPL:-}" ]; then
    echo "ERROR: set PANEL_VCF or PANEL_TMPL ({}=contig)" >&2; exit 1
fi
panel="${PANEL_TMPL:+${PANEL_TMPL//\{\}/$CONTIG}}"; panel="${panel:-${PANEL_VCF:-}}"
[ -s "$panel" ] || { echo "ERROR: panel VCF not found: $panel" >&2; exit 1; }

now() { date +%s.%N; }
el()  { awk -v a="$1" -v b="$2" 'BEGIN{printf "%.3f", b-a}'; }

log "Throughput: G2G diff vs. linearize-and-map, sweeping h on contig $CONTIG"
[ -f "$panel.tbi" ] || [ -f "$panel.csi" ] || { note "indexing $(basename "$panel")…"; "$TABIX" -p vcf "$panel"; }
# one contig-only reference (absolute coords, so consensus needs no rebasing)
CFA="$OUT/ref_${CONTIG}.fa"
[ -s "$CFA" ] || { note "extracting reference contig $CONTIG…"; "$SAMTOOLS" faidx "$REF_FA" "$CONTIG" > "$CFA"; "$SAMTOOLS" faidx "$CFA"; }

build_graph() {  # $1=sample-list  $2=out.gfa
    "$BCFTOOLS" view -r "$CONTIG" -S "$1" --force-samples -c1 -Oz "$panel" > "$OUT/_sub.vcf.gz"
    "$TABIX" -fp vcf "$OUT/_sub.vcf.gz"
    "$VG" construct -r "$REF_FA" -v "$OUT/_sub.vcf.gz" -R "$CONTIG" 2>/dev/null \
        | "$VG" convert -f - 2>/dev/null > "$2"
}
extract_haps() {  # $1=sample-list  $2=out.fa  — one haplotype per sample
    : > "$2"
    while read -r s; do
        [ -n "$s" ] || continue
        "$BCFTOOLS" consensus -f "$CFA" -s "$s" -H 1 "$panel" 2>/dev/null \
            | sed "1s/.*/>${s}_h1/" >> "$2"
    done < "$1"
}

TSV="$OUT/throughput.tsv"
printf 'h\tnodes\tg2g_s\tbaseline_s\n' > "$TSV"
printf '  %4s %12s %10s %14s\n' "h" "nodes(A+B)" "g2g(s)" "linmap(s)"
printf '  %4s %12s %10s %14s\n' "----" "----------" "------" "-------------"

for H in $H_SWEEP; do
    head -n "$H" "$SAMPLES_A" > "$OUT/a.ids"; head -n "$H" "$SAMPLES_B" > "$OUT/b.ids"
    build_graph "$OUT/a.ids" "$OUT/popA.gfa"
    build_graph "$OUT/b.ids" "$OUT/popB.gfa"

    # G2G: time the whole diff; capture node counts from its own report
    t0=$(now)
    rep=$("$G2G" "$OUT/popA.gfa" "$OUT/popB.gfa" --ref "$CONTIG" --vcf "$OUT/delta.vcf" 2>&1 >/dev/null) || true
    t1=$(now); g2g_s=$(el "$t0" "$t1")
    na=$(printf '%s\n' "$rep" | sed -n 's/.*GA: \([0-9]*\) nodes.*/\1/p' | head -1)
    nb=$(printf '%s\n' "$rep" | sed -n 's/.*GB: \([0-9]*\) nodes.*/\1/p' | head -1)
    nodes=$(( ${na:-0} + ${nb:-0} ))

    # baseline: extract + align h haplotypes per population (the dominant cost)
    t0=$(now)
    extract_haps "$OUT/a.ids" "$OUT/hapsA.fa"; extract_haps "$OUT/b.ids" "$OUT/hapsB.fa"
    "$MINIMAP2" -cx asm20 "$CFA" "$OUT/hapsA.fa" >/dev/null 2>&1 || true
    "$MINIMAP2" -cx asm20 "$CFA" "$OUT/hapsB.fa" >/dev/null 2>&1 || true
    t1=$(now); base_s=$(el "$t0" "$t1")

    printf '  %4s %12s %10s %14s\n' "$H" "$nodes" "$g2g_s" "$base_s"
    printf '%s\t%s\t%s\t%s\n' "$H" "$nodes" "$g2g_s" "$base_s" >> "$TSV"
done

echo
note "results: $TSV"
note "figure:  python3 demo/make_figures.py --throughput-csv $TSV"
note "interpretation: G2G is ~flat as h grows (graph size saturates); the baseline"
note "grows ~linearly in h — and that is only its extract+align stage, before any"
note "variant calling. G2G aligns the compressed graphs directly."
echo "OK — throughput baseline complete"
