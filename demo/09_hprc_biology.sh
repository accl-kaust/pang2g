#!/usr/bin/env bash
# Demo 6 — Q4 biological validation: are G2G delta alleles enriched for
#          population-differentiated structural variation?
#
# The claim under test: when G2G diffs two *population* pangenomes (here AFR vs
# EAS), the alleles it reports as present in one but not the other — the delta
# graph — should sit closer than chance to structural variants that are highly
# differentiated between those same populations (high Hudson F_ST in gnomAD-SV),
# and should also touch known pathogenic (ClinVar) sites. We test that with a
# permutation null (demo/hprc_enrich.py).
#
# Two modes:
#   --synthetic   End-to-end wiring + method check that runs ANYWHERE the g2g
#                 binary is built (no vg/bcftools/HPRC download needed). It plants
#                 a known population difference, runs the REAL g2g --ref diff, then
#                 synthesises a gnomAD-SV/ClinVar fixture and confirms the pipeline
#                 (a) recovers the planted enrichment and (b) reports NO enrichment
#                 for a shuffled-delta negative control. Produces the real TSV that
#                 figure fig_biology reads.
#   (default)     Real cluster pipeline: build AFR/EAS pangenomes with vg construct
#                 from a population-labelled panel, diff with g2g --ref, and enrich
#                 against real gnomAD-SV + ClinVar. Requires vg, bcftools, tabix and
#                 the data files below (point the script at them via env vars).
#
# Real-mode inputs (env vars; any left unset are downloaded or, if impossible,
# cause a clear error):
#   REF_FA        GRCh38 fasta (+ .fai)          CHROMS        e.g. "chr20 chr21 chr22"
#   PANEL_VCF     population-labelled phased VCF  SAMPLES_A/_B  sample-id files (one/line)
#   GNOMAD_SV_VCF gnomAD-SV sites VCF (AFR_AF,EAS_AF)            CLINVAR_VCF   ClinVar VCF
#   CHROM_SIZES   GRCh38 .fai or .chrom.sizes
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
need_g2g
OUT="$WORK/biology"; mkdir -p "$OUT"
ENRICH="$DEMO_DIR/hprc_enrich.py"
WINDOW="${WINDOW:-400}"; PERMS="${PERMS:-1000}"
FST="${FST_THRESHOLD:-0.2}"; AFR="${AFR_FIELD:-AFR_AF}"; EAS="${EAS_FIELD:-EAS_AF}"

# ── helper: assert a result row passes (fold and p within expected band) ───────
check() {  # $1=tsv $2=row-regex $3=label $4=min_fold $5=max_fold $6=max_p
    local row fold p
    row=$(awk -F'\t' -v re="$2" 'NR>1 && $1 ~ re {print; exit}' "$1")
    [ -z "$row" ] && { echo "  [FAIL] $3: no '$2' row in $1"; return 1; }
    fold=$(echo "$row" | cut -f4); p=$(echo "$row" | cut -f5)
    awk -v f="$fold" -v p="$p" -v lo="$4" -v hi="$5" -v mp="$6" -v lab="$3" 'BEGIN{
        ok = (f+0 >= lo && f+0 <= hi && p+0 <= mp)
        printf "  [%s] %-34s fold=%.2fx  p=%s  (want %g..%g, p<=%g)\n",
               (ok?"PASS":"FAIL"), lab, f, p, lo, hi, mp
        exit !ok }'
}

# ══════════════════════════════════════════════════════════════════════════════
#  SYNTHETIC MODE — runs with only g2g + python stdlib (validated in CI/laptop)
# ══════════════════════════════════════════════════════════════════════════════
run_synthetic() {
    local S="$OUT/synthetic"; mkdir -p "$S"
    local LEN="${SYN_LEN:-10000000}" CHROM="chr1"
    log "Synthetic Q4 — plant a population difference, diff it, test enrichment"

    note "1. build two population pangenomes that differ at known coordinates"
    python3 "$DEMO_DIR/make_synthetic_pair.py" \
        --length "$LEN" --n-shared 300 --n-a 140 --n-b 160 \
        --chop-a 32 --chop-b 25 --ref-name "$CHROM" --outdir "$S" >/dev/null
    note "   $(grep -cP '\t[AB]$' "$S/truth.tsv") population-specific sites planted on $CHROM"

    note "2. g2g --ref $CHROM  (recover the delta alleles at real coordinates)"
    "$G2G" "$S/demoA.gfa" "$S/demoB.gfa" --ref "$CHROM" --vcf "$S/delta.vcf" 2>&1 \
        | grep -o "a_only=[0-9]*  b_only=[0-9]*" | sed 's/^/     /'

    note "3. synthesise a gnomAD-SV + ClinVar fixture around the recovered deltas"
    # High-F_ST SVs near 60% of deltas (true signal) PLUS an equal number of
    # high-F_ST decoys at random loci (so a delta can MISS one — the test is not
    # rigged to overlap), PLUS 4x neutral low-F_ST SVs (exercise the F_ST filter).
    python3 - "$S/delta.vcf" "$CHROM" "$LEN" "$S" 7 "$WINDOW" <<'PY'
import sys, random
delta, chrom, length, outdir, seed, window = sys.argv[1:7]
length, seed, window = int(length), int(seed), int(window)
rng = random.Random(seed)
deltas = [(c.split('\t')[0], int(c.split('\t')[1]))
          for c in open(delta) if not c.startswith('#') and len(c.split('\t')) >= 8]
hit = rng.sample(deltas, max(1, int(0.6 * len(deltas))))
recs = []                                              # (pos, AFR_AF, EAS_AF)
for _, p in hit:                                       # high-F_ST, near a delta
    recs.append((max(1, p + rng.randint(-window // 3, window // 3)), 0.04, 0.88))
for _ in range(len(hit)):                              # high-F_ST decoys, random
    recs.append((rng.randint(1, length), 0.05, 0.86))
for _ in range(4 * len(deltas)):                       # neutral, low F_ST
    a = round(rng.uniform(0.20, 0.45), 3)
    recs.append((rng.randint(1, length), a, round(max(0.01, a + rng.uniform(-0.05, 0.05)), 3)))
recs.sort()
with open(f"{outdir}/gnomad_sv.vcf", "w") as f:
    f.write(f"##fileformat=VCFv4.2\n##contig=<ID={chrom}>\n")
    f.write('##INFO=<ID=AFR_AF,Number=1,Type=Float,Description="AFR allele freq">\n')
    f.write('##INFO=<ID=EAS_AF,Number=1,Type=Float,Description="EAS allele freq">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for i, (p, a, b) in enumerate(recs):
        f.write(f"{chrom}\t{p}\tsv{i}\tN\t<INS>\t.\tPASS\tSVTYPE=INS;AFR_AF={a};EAS_AF={b}\n")
cv = [(max(1, p + rng.randint(-window // 3, window // 3)), "Pathogenic")
      for _, p in rng.sample(deltas, max(1, int(0.25 * len(deltas))))]
cv += [(rng.randint(1, length), "Benign") for _ in range(3 * len(deltas))]
cv.sort()
with open(f"{outdir}/clinvar.vcf", "w") as f:
    f.write(f"##fileformat=VCFv4.2\n##contig=<ID={chrom}>\n")
    f.write('##INFO=<ID=CLNSIG,Number=.,Type=String,Description="clinical significance">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for i, (p, sig) in enumerate(cv):
        f.write(f"{chrom}\t{p}\tcv{i}\tA\tT\t.\tPASS\tCLNSIG={sig}\n")
with open(f"{outdir}/chr.sizes", "w") as f:
    f.write(f"{chrom}\t{length}\n")
print(f"     gnomAD-SV: {len(recs)} records ({2*len(hit)} high-F_ST), "
      f"ClinVar: {len(cv)} ({len(cv)-3*len(deltas)} pathogenic)")
PY

    note "4. enrichment test: delta vs high-F_ST gnomAD-SV + ClinVar (permutation null)"
    python3 "$ENRICH" --delta "$S/delta.vcf" --gnomad-sv "$S/gnomad_sv.vcf" \
        --sizes "$S/chr.sizes" --clinvar "$S/clinvar.vcf" \
        --afr-field AFR_AF --eas-field EAS_AF --fst-threshold "$FST" \
        --window "$WINDOW" --perms "$PERMS" --seed 11 \
        --out "$OUT/biology_results.tsv"

    note "5. negative control: shuffle the deltas to random loci → enrichment must vanish"
    python3 - "$S/delta.vcf" "$LEN" "$S/delta_shuf.vcf" <<'PY'
import sys, random
src, length, dst = sys.argv[1], int(sys.argv[2]), sys.argv[3]
rng = random.Random(99); out = []
for L in open(src):
    if L.startswith('#'): out.append(L); continue
    c = L.rstrip('\n').split('\t')
    if len(c) >= 8: c[1] = str(rng.randint(1, length)); out.append('\t'.join(c) + '\n')
open(dst, 'w').writelines(out)
PY
    python3 "$ENRICH" --delta "$S/delta_shuf.vcf" --gnomad-sv "$S/gnomad_sv.vcf" \
        --sizes "$S/chr.sizes" --afr-field AFR_AF --eas-field EAS_AF \
        --fst-threshold "$FST" --window "$WINDOW" --perms "$PERMS" --seed 13 \
        --out "$S/control.tsv" >/dev/null

    echo
    log "Synthetic validation summary"
    local rc=0
    check "$OUT/biology_results.tsv" "F_ST"    "delta vs high-F_ST (signal)"  1.5 100000 0.05 || rc=1
    check "$OUT/biology_results.tsv" "ClinVar" "delta vs ClinVar pathogenic"  1.3 100000 0.05 || rc=1
    check "$S/control.tsv"           "F_ST"    "shuffled control (null)"      0.3 2.5    1.00 || rc=1
    echo
    note "results TSV: $OUT/biology_results.tsv"
    note "figure:  python3 demo/make_figures.py --biology-csv $OUT/biology_results.tsv"
    [ "$rc" = 0 ] && echo "OK — synthetic Q4 pipeline validated (signal enriched, control null)" \
                  || { echo "FAILED — see rows above"; return 1; }
}

# ══════════════════════════════════════════════════════════════════════════════
#  REAL MODE — population pangenomes from a labelled panel (needs vg + bcftools)
# ══════════════════════════════════════════════════════════════════════════════
run_real() {
    for t in vg bcftools tabix; do
        have "$t" || { echo "ERROR: '$t' not found. Install it, or run with --synthetic." >&2; exit 1; }
    done
    : "${REF_FA:?set REF_FA to a GRCh38 fasta (with .fai)}"
    : "${PANEL_VCF:?set PANEL_VCF to a population-labelled phased VCF}"
    : "${SAMPLES_A:?set SAMPLES_A to a file of AFR sample ids (one per line)}"
    : "${SAMPLES_B:?set SAMPLES_B to a file of EAS sample ids (one per line)}"
    : "${GNOMAD_SV_VCF:?set GNOMAD_SV_VCF to a gnomAD-SV sites VCF (AFR_AF,EAS_AF)}"
    : "${CHROM_SIZES:?set CHROM_SIZES to a GRCh38 .fai or .chrom.sizes}"
    local CHROMS="${CHROMS:-chr20 chr21 chr22}"
    local R="$OUT/real"; mkdir -p "$R"

    log "Real Q4 — build AFR/EAS pangenomes, diff with g2g --ref, enrich"
    : > "$R/delta.vcf.body"
    for chrom in $CHROMS; do
        note "[$chrom] subset panel by population → vg construct → GFA"
        for grp in A B; do
            local samp; samp=$([ "$grp" = A ] && echo "$SAMPLES_A" || echo "$SAMPLES_B")
            bcftools view -r "$chrom" -S "$samp" --force-samples -c1 -Oz \
                "$PANEL_VCF" > "$R/pop${grp}_${chrom}.vcf.gz"
            tabix -fp vcf "$R/pop${grp}_${chrom}.vcf.gz"
            vg construct -r "$REF_FA" -v "$R/pop${grp}_${chrom}.vcf.gz" -R "$chrom" 2>/dev/null \
                | vg convert -f - 2>/dev/null > "$R/pop${grp}_${chrom}.gfa"
        done
        note "[$chrom] g2g --ref $chrom"
        "$G2G" "$R/popA_${chrom}.gfa" "$R/popB_${chrom}.gfa" --ref "$chrom" \
            --vcf "$R/delta_${chrom}.vcf"
        grep -v '^#' "$R/delta_${chrom}.vcf" >> "$R/delta.vcf.body" || true
    done
    { echo "##fileformat=G2GdeltaVCF-0.2";
      echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
      cat "$R/delta.vcf.body"; } > "$R/delta.vcf"
    note "total delta alleles: $(grep -vc '^#' "$R/delta.vcf")"

    log "Enrichment test against real gnomAD-SV + ClinVar"
    python3 "$ENRICH" --delta "$R/delta.vcf" --gnomad-sv "$GNOMAD_SV_VCF" \
        --sizes "$CHROM_SIZES" ${CLINVAR_VCF:+--clinvar "$CLINVAR_VCF"} \
        --afr-field "$AFR" --eas-field "$EAS" --fst-threshold "$FST" \
        --window "$WINDOW" --perms "$PERMS" --out "$OUT/biology_results.tsv"
    echo; note "results TSV: $OUT/biology_results.tsv"
    note "figure:  python3 demo/make_figures.py --biology-csv $OUT/biology_results.tsv"
    echo "OK — real Q4 pipeline complete"
}

case "${1:-}" in
    --synthetic) run_synthetic ;;
    --self-test) python3 "$ENRICH" --self-test ;;
    "")          run_real ;;
    *) echo "usage: $0 [--synthetic | --self-test]   (no arg = real cluster pipeline)"; exit 2 ;;
esac
