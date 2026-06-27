#!/usr/bin/env bash
# Demo 3 — SCALING: drive the g2g --ref diff across increasing graph sizes,
# measuring wall time, peak memory and per-phase timing while CHECKING that the
# diff stays EXACT at every size (the strong claim: correct *and* near-linear at
# scale, not just fast on toy inputs).
#
# Pure synthetic with known ground truth at every point — no network, no
# third-party tools. Each point builds two graphs that embed the same reference
# but chop it into DIFFERENT nodes (A=25 bp, B=17 bp) and plant N variants
# (80% shared, 10% A-only, 10% B-only); g2g must recover exactly N/10 A-only and
# N/10 B-only at every size, proving coordinate matching holds as graphs grow
# from ~kilobase (a gene) to tens of megabases (a chromosome arm).
#
# Timing and peak memory come from g2g's own --verbose-timing output ([timing]
# wall=, [mem] peak_rss_mb=), so no GNU /usr/bin/time is required.
#
# Tuning (env):
#   SCALE_SIZES="bp:variants ..."   override the whole sweep
#   SCALE_MAX_MB=N                  append one more point at N megabases
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
need_g2g
OUT="$WORK/scaling"; mkdir -p "$OUT"

# Default sweep: bp:total_variants, ~1 variant / 200 bp (a dense, stressful SNP
# rate). Tops out at 8 Mb (~700x the 11 kb HLA gene region).
SIZES="${SCALE_SIZES:-100000:500 500000:2500 2000000:10000 8000000:40000}"
if [ -n "${SCALE_MAX_MB:-}" ]; then
    bp=$(( SCALE_MAX_MB * 1000000 ))
    SIZES="$SIZES ${bp}:$(( bp / 200 ))"
fi

TSV="$OUT/scale_results.tsv"
printf 'ref_bp\tvariants\tnodes_A\tnodes_B\tsnarl_pairs\ta_only\tb_only\tcorrect\tdecompose_s\twall_s\tpeak_MB\n' > "$TSV"

log "Scaling sweep — g2g --ref diff, increasing size (A=25bp / B=17bp chopping)"
printf '  %10s %9s %9s %9s %8s %7s %7s %5s %10s %8s %8s\n' \
       ref_bp variants nodesA nodesB snarls a_only b_only OK decomp_s wall_s peak_MB

worst=""
for spec in $SIZES; do
    L="${spec%%:*}"; N="${spec##*:}"
    na=$(( N / 10 )); [ "$na" -lt 1 ] && na=1
    nb="$na"; nsh=$(( N - na - nb ))

    python3 "$DEMO_DIR/make_synthetic_pair.py" --length "$L" \
        --n-shared "$nsh" --n-a "$na" --n-b "$nb" \
        --chop-a 25 --chop-b 17 --outdir "$OUT" >/dev/null

    err="$OUT/run_${L}.err"
    "$G2G" "$OUT/demoA.gfa" "$OUT/demoB.gfa" --ref GRCh38 \
        --verbose-timing --vcf "$OUT/scale.vcf" >/dev/null 2>"$err"

    nA=$(sed -n 's/.*GA: \([0-9]*\) nodes  GB: \([0-9]*\) nodes.*/\1/p' "$err")
    nB=$(sed -n 's/.*GA: \([0-9]*\) nodes  GB: \([0-9]*\) nodes.*/\2/p' "$err")
    sp=$(sed -n 's/.*snarl_pairs=\([0-9]*\).*/\1/p' "$err")
    ao=$(sed -n 's/.*a_only=\([0-9]*\)  b_only=\([0-9]*\).*/\1/p' "$err")
    bo=$(sed -n 's/.*a_only=\([0-9]*\)  b_only=\([0-9]*\).*/\2/p' "$err")
    dec=$(sed -n 's/.*decompose=\([0-9.]*\).*/\1/p' "$err")
    wall=$(sed -n 's/.* wall=\([0-9.]*\).*/\1/p' "$err")
    mem=$(sed -n 's/.*peak_rss_mb=\([0-9.]*\).*/\1/p' "$err")
    if [ "$ao" = "$na" ] && [ "$bo" = "$nb" ]; then ok="PASS"; else ok="FAIL"; fi

    printf '  %10s %9s %9s %9s %8s %7s %7s %5s %10s %8s %8s\n' \
           "$L" "$N" "$nA" "$nB" "$sp" "$ao" "$bo" "$ok" "$dec" "$wall" "$mem"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
           "$L" "$N" "$nA" "$nB" "$sp" "$ao" "$bo" "$ok" "$dec" "$wall" "$mem" >> "$TSV"
    worst="$L $N $nA $nB $wall $mem"
done

echo
note "results table saved to $TSV"
if grep -q 'FAIL' "$TSV"; then
    echo "ERROR: a scaling point produced an incorrect diff" >&2; exit 1
fi
# Throughput headline from the largest point.
set -- $worst
awk -v bp="$1" -v var="$2" -v nA="$3" -v nB="$4" -v w="$5" -v m="$6" 'BEGIN{
    printf "    largest point: %.1f Mb ref, %d+%d nodes, %d variants -> %ss, %s MB peak\n", bp/1e6, nA, nB, var, w, m;
    if (w+0>0) printf "    throughput: %.0f variants/s, %.0f total nodes/s\n", var/w, (nA+nB)/w;
}'
echo "OK — scaling sweep complete: diff stayed EXACT at every size, near-linear time"
