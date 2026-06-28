#!/usr/bin/env bash
# Demo 4 — parallel scaling of the work-stealing aligner.
#
# Fixes one graph pair and sweeps the thread count, reporting end-to-end AND
# DP-phase speedup (the Amdahl distinction) plus the scheduler's [runstats]
# (tasks / steals / imbalance). Writes a TSV that make_figures.py turns into
# fig_scaling.
#
# What is measured, and why these choices:
#   * Baseline is --threads 1, NOT --sequential. Both use the SAME algorithm
#     (anchor-based snarl correspondence + per-snarl DP); the only difference is
#     thread count, so the speedup is pure parallelism. (--sequential is a
#     DIFFERENT, greedy O(snarls^2) pairing — it over-pairs and is far slower;
#     comparing against it conflates an algorithmic change with parallelism, as
#     an early "172x" measurement did.)
#   * The pair carries SV-like snarls (long, multi-allele bubbles via --alt-len),
#     because G2G parallelises the per-snarl matching DP. With SNP-only bubbles
#     that DP is trivial and dominated by the (sequential) decompose/I-O, so no
#     workload reaches the parallel engine. Long alleles make the DP the bulk of
#     the work, which is the regime real structural-variant regions live in.
#   * Equal chopping (CHOP for both) so the anchor correspondence has matching
#     backbone nodes.
#
# Env:
#   SIZE_MB  reference size of the benchmark pair        (default 8)
#   CHOP     node size for BOTH graphs                   (default 25)
#   MAX_K    max alleles per snarl                       (default 40)
#   CFRAC    fraction of shared sites made multi-allele  (default 0.5)
#   ALT_LEN  max length of alt alleles (SV-like = long)  (default 200)
#   THREADS  thread counts to sweep   (default "1 2 4 8 16 32 64 128")
#   REPS     repetitions; the min wall is kept           (default 3)
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
need_g2g
OUT="$WORK/parallel"; mkdir -p "$OUT"
SIZE_MB="${SIZE_MB:-8}"; CHOP="${CHOP:-25}"; MAX_K="${MAX_K:-40}"
CFRAC="${CFRAC:-0.5}"; ALT_LEN="${ALT_LEN:-200}"; REPS="${REPS:-3}"
NPROC="$(nproc)"; THREADS="${THREADS:-1 2 4 8 16 32 64 128}"

bp=$(( SIZE_MB * 1000000 )); N=$(( bp / 200 )); na=$(( N / 10 )); nb="$na"; nsh=$(( N - na - nb ))
log "Parallel scaling — ${SIZE_MB} Mb pair, SV-like snarls (max_k=${MAX_K}, alt_len=${ALT_LEN}); baseline = 1 thread"
python3 "$DEMO_DIR/make_synthetic_pair.py" --length "$bp" \
    --n-shared "$nsh" --n-a "$na" --n-b "$nb" --chop-a "$CHOP" --chop-b "$CHOP" \
    --max-k "$MAX_K" --complex-frac "$CFRAC" --alt-len "$ALT_LEN" --outdir "$OUT" >/dev/null
note "GA/GB segments: $(grep -c '^S' "$OUT/demoA.gfa") / $(grep -c '^S' "$OUT/demoB.gfa")"

run() {  # $1 = threads; echoes "wall dp tasks steals imbalance pairs"
    local t="$1" err w best_w="" best=""
    for _ in $(seq 1 "$REPS"); do
        err=$("$G2G" "$OUT/demoA.gfa" "$OUT/demoB.gfa" --threads "$t" --verbose-timing 2>&1 >/dev/null) || true
        w=$(printf '%s\n' "$err" | sed -n 's/.* wall=\([0-9.]*\).*/\1/p')
        [ -z "$w" ] && continue
        if [ -z "$best_w" ] || awk -v a="$w" -v b="$best_w" 'BEGIN{exit !(a<b)}'; then best_w="$w"; best="$err"; fi
    done
    local dp tk se im pp
    dp=$(printf '%s\n' "$best" | sed -n 's/.* dp=\([0-9.]*\).*/\1/p')
    tk=$(printf '%s\n' "$best" | sed -n 's/.*tasks=\([0-9]*\).*/\1/p')
    se=$(printf '%s\n' "$best" | sed -n 's/.*steals=\([0-9]*\).*/\1/p')
    im=$(printf '%s\n' "$best" | sed -n 's/.*imbalance=\([0-9.]*\).*/\1/p')
    pp=$(printf '%s\n' "$best" | sed -n 's/.*snarl_pairs=\([0-9]*\).*/\1/p')
    echo "${best_w:-NA} ${dp:-NA} ${tk:-NA} ${se:-NA} ${im:-NA} ${pp:-NA}"
}
spd() { awk -v a="$1" -v b="$2" 'BEGIN{ if (b>0) printf "%.2f", a/b; else printf "NA" }'; }

TSV="$OUT/parallel_scaling.tsv"
printf 'threads\twall_s\tdp_s\ttasks\tsteals\timbalance\tpairs\tspeedup_e2e\tspeedup_dp\tconsistent\n' > "$TSV"
printf '  %8s %9s %9s %10s %9s %10s %7s %8s %7s\n' threads wall dp tasks steals imbal pairs e2e dp

B1W=""; B1D=""; BP=""; bad=0
for t in $THREADS; do
    if [ "$t" -gt "$NPROC" ]; then continue; fi
    read -r W D TK SE IM P < <(run "$t")
    if [ -z "$B1W" ]; then B1W="$W"; B1D="$D"; BP="$P"; fi   # baseline = first (1 thread)
    e=$(spd "$B1W" "$W"); dps=$(spd "$B1D" "$D")
    cons="ok"; [ "${P:-x}" != "$BP" ] && { cons="PAIRS≠"; bad=1; }   # all thread counts must agree
    printf '  %8s %9s %9s %10s %9s %10s %7s %8s %7s\n' "$t" "$W" "${D:-NA}" "${TK:-NA}" "${SE:-NA}" "${IM:-NA}" "${P:-NA}" "$e" "$dps"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$t" "$W" "${D:-NA}" "${TK:-NA}" "${SE:-NA}" "${IM:-NA}" "${P:-NA}" "$e" "$dps" "$cons" >> "$TSV"
done

echo
note "results: $TSV   (baseline = 1 thread; speedup is pure parallelism)"
note "plot:    python3 make_figures.py --scaling-csv $TSV"
if [ "$bad" = 1 ]; then
    echo "WARNING: thread counts disagree on snarl-pair count — a parallel bug, not a speedup." >&2
else
    echo "OK — parallel scaling sweep complete (deterministic across thread counts)"
fi
