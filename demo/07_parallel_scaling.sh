#!/usr/bin/env bash
# Demo 4 — parallel scaling of the work-stealing aligner.
#
# Fixes one graph pair (heavy-tailed multi-allele snarls, so the per-snarl DP
# carries real, skewed work) and sweeps the thread count, reporting end-to-end
# AND DP-phase speedup (the Amdahl distinction the paper makes) plus the
# scheduler's [runstats] (tasks / steals / imbalance). Writes a TSV that
# make_figures.py turns into fig_scaling.
#
# This benchmarks the snarl-pairing ALIGNER (g2g A B, no --ref) — a different
# engine from the single-threaded reference diff of Demos 3/3b. The pair is
# built with EQUAL chopping (CHOP for both A and B) because the aligner pairs
# snarls via shared backbone anchors; unequal chopping leaves no matching
# backbone nodes and the aligner finds nothing.
#
# IMPORTANT: every run is checked against the sequential aligner. If the
# parallel result disagrees (snarl-pair count or total score), the row is
# flagged DIVERGES and the speedup is meaningless — fix the aligner first.
#
# Env:
#   SIZE_MB  reference size of the benchmark pair        (default 32)
#   CHOP     node size for BOTH graphs                   (default 25)
#   MAX_K    max alleles per snarl (heavy tail)          (default 40)
#   CFRAC    fraction of shared sites made multi-allele  (default 0.5)
#   THREADS  thread counts to sweep   (default "1 2 4 8 16 32 64 128")
#   REPS     repetitions; the min wall is kept           (default 3)
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
need_g2g
OUT="$WORK/parallel"; mkdir -p "$OUT"
SIZE_MB="${SIZE_MB:-32}"; CHOP="${CHOP:-25}"; MAX_K="${MAX_K:-40}"; CFRAC="${CFRAC:-0.5}"
REPS="${REPS:-3}"; NPROC="$(nproc)"; THREADS="${THREADS:-1 2 4 8 16 32 64 128}"

bp=$(( SIZE_MB * 1000000 )); N=$(( bp / 200 )); na=$(( N / 10 )); nb="$na"; nsh=$(( N - na - nb ))
log "Parallel scaling — ${SIZE_MB} Mb pair, equal chop=${CHOP}, heavy-tailed snarls (max_k=${MAX_K}); sweep threads"
python3 "$DEMO_DIR/make_synthetic_pair.py" --length "$bp" \
    --n-shared "$nsh" --n-a "$na" --n-b "$nb" --chop-a "$CHOP" --chop-b "$CHOP" \
    --max-k "$MAX_K" --complex-frac "$CFRAC" --outdir "$OUT" >/dev/null
note "GA/GB segments: $(grep -c '^S' "$OUT/demoA.gfa") / $(grep -c '^S' "$OUT/demoB.gfa")"

# Run the aligner; echo "wall dp tasks steals imbalance pairs score" (NA if absent).
run() {  # $1 = "seq" | <threads>
    local mode="$1" err w best_w="" best=""
    for _ in $(seq 1 "$REPS"); do
        if [ "$mode" = seq ]; then
            err=$("$G2G" "$OUT/demoA.gfa" "$OUT/demoB.gfa" --sequential --verbose-timing 2>&1 >/dev/null) || true
        else
            err=$("$G2G" "$OUT/demoA.gfa" "$OUT/demoB.gfa" --threads "$mode" --verbose-timing 2>&1 >/dev/null) || true
        fi
        w=$(printf '%s\n' "$err" | sed -n 's/.* wall=\([0-9.]*\).*/\1/p')
        [ -z "$w" ] && continue
        if [ -z "$best_w" ] || awk -v a="$w" -v b="$best_w" 'BEGIN{exit !(a<b)}'; then
            best_w="$w"; best="$err"
        fi
    done
    # Capture into vars and default missing fields to NA, so the consumer's
    # `read` always sees exactly 7 columns (the sequential run has no [runstats]).
    local dp tk se im pp sc
    dp=$(printf '%s\n' "$best" | sed -n 's/.* dp=\([0-9.]*\).*/\1/p')
    tk=$(printf '%s\n' "$best" | sed -n 's/.*tasks=\([0-9]*\).*/\1/p')
    se=$(printf '%s\n' "$best" | sed -n 's/.*steals=\([0-9]*\).*/\1/p')
    im=$(printf '%s\n' "$best" | sed -n 's/.*imbalance=\([0-9.]*\).*/\1/p')
    pp=$(printf '%s\n' "$best" | sed -n 's/.*snarl_pairs=\([0-9]*\).*/\1/p')
    sc=$(printf '%s\n' "$best" | sed -n 's/.*total_score=\([0-9.]*\).*/\1/p')
    echo "${best_w:-NA} ${dp:-NA} ${tk:-NA} ${se:-NA} ${im:-NA} ${pp:-NA} ${sc:-NA}"
}
spd() { awk -v a="$1" -v b="$2" 'BEGIN{ if (b>0) printf "%.2f", a/b; else printf "NA" }'; }

TSV="$OUT/parallel_scaling.tsv"
printf 'threads\twall_s\tdp_s\ttasks\tsteals\timbalance\tpairs\tscore\tconsistent\tspeedup_e2e\tspeedup_dp\n' > "$TSV"
printf '  %7s %9s %9s %10s %9s %10s %7s %6s %8s %7s\n' threads wall dp tasks steals imbal pairs cons e2e dp
read -r SW SD ST SS SI SP SC < <(run seq)
printf '  %7s %9s %9s %10s %9s %10s %7s %6s %8s %7s\n' seq "$SW" "${SD:-NA}" "${ST:-NA}" "${SS:-NA}" "${SI:-NA}" "${SP:-NA}" "ref" "1.00" "1.00"
printf 'seq\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tref\t1.00\t1.00\n' "$SW" "${SD:-NA}" "${ST:-NA}" "${SS:-NA}" "${SI:-NA}" "${SP:-NA}" "${SC:-NA}" >> "$TSV"

diverged=0
for t in $THREADS; do
    if [ "$t" -gt "$NPROC" ]; then continue; fi
    read -r W D TK SE IM P C < <(run "$t")
    # consistency vs sequential: same pair count AND score within 0.1%
    cons="ok"
    if [ "${P:-NA}" != "$SP" ] || awk -v a="${C:-0}" -v b="${SC:-0}" 'BEGIN{exit !(b>0 && (a-b<-0.001*b || a-b>0.001*b))}'; then
        cons="DIVERGES"; diverged=1
    fi
    e=$(spd "$SW" "$W"); dps=$(spd "$SD" "$D")
    printf '  %7s %9s %9s %10s %9s %10s %7s %6s %8s %7s\n' "$t" "$W" "${D:-NA}" "${TK:-NA}" "${SE:-NA}" "${IM:-NA}" "${P:-NA}" "$cons" "$e" "$dps"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$t" "$W" "${D:-NA}" "${TK:-NA}" "${SE:-NA}" "${IM:-NA}" "${P:-NA}" "${C:-NA}" "$cons" "$e" "$dps" >> "$TSV"
done

echo
note "results: $TSV"
if [ "$diverged" = 1 ]; then
    echo "WARNING: the parallel aligner DISAGREES with the sequential aligner on this"  >&2
    echo "         graph (snarl-pair count / total score differ). The speedups above"   >&2
    echo "         are therefore NOT a like-for-like comparison — reconcile the parallel" >&2
    echo "         correspondence with the sequential aligner before reporting them."     >&2
else
    note "plot: python3 make_figures.py --scaling-csv $TSV   (fills fig_scaling)"
    echo "OK — parallel scaling sweep complete (parallel == sequential at every thread count)"
fi
