#!/usr/bin/env bash
# Demo 3b - push the scaling to the machine's limit.
#
# A geometric sweep (doubling the reference size each step) that climbs until
# the NEXT point would exceed the node's RAM or a wall-clock budget, then stops
# cleanly and reports the ceiling. It projects each point's peak memory and wall
# time from the previous one and refuses to *attempt* a point it expects to OOM,
# so an unattended run won't get killed and waste a cluster allocation. The diff
# is checked for exact correctness at every size.
#
# Metrics are g2g's own (--verbose-timing: wall=, peak_rss_mb=), so no external
# profiler is needed. GFA files are deleted after each point to bound disk.
#
# Env knobs (all optional):
#   BUDGET_SEC=N    wall-clock budget in seconds (default 32400 = 9 h; keep it
#                   comfortably under your node's hard limit)
#   MEM_CEIL_MB=N   peak-RSS ceiling (default 80% of MemTotal)
#   START_MB=N      first reference size, megabases (default 8)
#   MAX_MB=N        hard cap on reference size (default 8192 = 8 Gb)
#   DENSITY=N       bp per planted variant (default 200)
#
# Recommended invocation (captures a log you can commit as evidence):
#   bash demo/06_scaling_limit.sh 2>&1 | tee demo/logs/06_scaling_limit.log
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
need_g2g
OUT="$WORK/scaling_limit"; mkdir -p "$OUT"

mt_kb=$(awk '/MemTotal/{print $2}' /proc/meminfo)
BUDGET_SEC="${BUDGET_SEC:-32400}"
MEM_CEIL_MB="${MEM_CEIL_MB:-$(awk -v k="$mt_kb" 'BEGIN{printf "%d", k/1024*0.80}')}"
START_MB="${START_MB:-8}"
MAX_MB="${MAX_MB:-8192}"
DENSITY="${DENSITY:-200}"
CHOP_A=25; CHOP_B=17

log "Scaling to the limit - doubling until RAM or budget is hit"
note "MemTotal $(awk -v k="$mt_kb" 'BEGIN{printf "%.1f", k/1048576}') GB  |  mem ceiling ${MEM_CEIL_MB} MB  |  budget ${BUDGET_SEC}s  |  density 1 var / ${DENSITY} bp"

TSV="$OUT/scale_limit_results.tsv"
printf 'ref_mb\tref_bp\tnodes_A\tnodes_B\tsnarls\ta_only\tb_only\tcorrect\tdecompose_s\twall_s\tpeak_MB\n' > "$TSV"
printf '  %7s %11s %11s %9s %8s %7s %5s %10s %9s %9s\n' \
       ref_mb nodesA nodesB snarls a_only OK decomp_s wall_s peak_MB elapsed_s

start=$(date +%s)
mb="$START_MB"
last_nodes=0; last_mem=0; last_wall=0; last_mb=0
ceiling=""; stop_reason=""

while :; do
    now=$(date +%s); elapsed=$(( now - start ))
    if [ "$mb" -gt "$MAX_MB" ]; then stop_reason="reached MAX_MB=${MAX_MB}"; break; fi

    bp=$(( mb * 1000000 ))
    N=$(( bp / DENSITY )); na=$(( N / 10 )); [ "$na" -lt 1 ] && na=1
    nb="$na"; nsh=$(( N - na - nb ))

    # -- projection guards (skip on the first, known-safe point) --------------
    if [ "$last_nodes" -gt 0 ]; then
        ratio=$(awk -v a="$mb" -v b="$last_mb" 'BEGIN{print a/b}')
        proj_mem=$(awk -v m="$last_mem" -v r="$ratio" 'BEGIN{printf "%.0f", m*r}')
        proj_wall=$(awk -v w="$last_wall" -v r="$ratio" 'BEGIN{printf "%.1f", w*(r^1.15)}')
        rem=$(( BUDGET_SEC - elapsed ))
        if awk -v p="$proj_mem" -v c="$MEM_CEIL_MB" 'BEGIN{exit !(p>c)}'; then
            stop_reason="next point ~${proj_mem} MB would exceed mem ceiling ${MEM_CEIL_MB} MB"; break; fi
        if awk -v p="$proj_wall" -v r="$rem" 'BEGIN{exit !(p>r)}'; then
            stop_reason="next point ~${proj_wall}s would exceed remaining budget ${rem}s"; break; fi
    fi

    # -- disk guard (GFA for both graphs, ~65 B/node) -------------------------
    est_mb=$(awk -v bp="$bp" -v ca="$CHOP_A" -v cb="$CHOP_B" 'BEGIN{printf "%.0f",(bp/ca+bp/cb)*65/1048576}')
    free_mb=$(df -Pm "$OUT" | awk 'NR==2{print $4}')
    if [ "$free_mb" -lt $(( est_mb * 2 )) ]; then
        stop_reason="insufficient disk for ~${est_mb} MB of GFA (free ${free_mb} MB)"; break; fi

    # -- generate -------------------------------------------------------------
    if ! python3 "$DEMO_DIR/make_synthetic_pair.py" --length "$bp" \
            --n-shared "$nsh" --n-a "$na" --n-b "$nb" \
            --chop-a "$CHOP_A" --chop-b "$CHOP_B" --outdir "$OUT" >/dev/null 2>"$OUT/gen.err"; then
        stop_reason="generator failed at ${mb} Mb ($(tail -1 "$OUT/gen.err"))"; break; fi

    # -- diff (robust: a kill/OOM here must not abort the script) -------------
    err="$OUT/run.err"
    if ! "$G2G" "$OUT/demoA.gfa" "$OUT/demoB.gfa" --ref GRCh38 \
            --verbose-timing --vcf "$OUT/scale.vcf" >/dev/null 2>"$err"; then
        rc=$?
        rm -f "$OUT"/demo*.gfa
        stop_reason="g2g exited ${rc} at ${mb} Mb (likely OOM-killed) - ceiling is the previous point"; break; fi

    nA=$(sed -n 's/.*GA: \([0-9]*\) nodes  GB: \([0-9]*\) nodes.*/\1/p' "$err")
    nB=$(sed -n 's/.*GA: \([0-9]*\) nodes  GB: \([0-9]*\) nodes.*/\2/p' "$err")
    sp=$(sed -n 's/.*snarl_pairs=\([0-9]*\).*/\1/p' "$err")
    ao=$(sed -n 's/.*a_only=\([0-9]*\)  b_only=\([0-9]*\).*/\1/p' "$err")
    bo=$(sed -n 's/.*a_only=\([0-9]*\)  b_only=\([0-9]*\).*/\2/p' "$err")
    dec=$(sed -n 's/.*decompose=\([0-9.]*\).*/\1/p' "$err")
    wall=$(sed -n 's/.* wall=\([0-9.]*\).*/\1/p' "$err")
    mem=$(sed -n 's/.*peak_rss_mb=\([0-9.]*\).*/\1/p' "$err")
    rm -f "$OUT"/demo*.gfa

    if [ -z "$wall" ] || [ -z "$mem" ] || [ -z "$nA" ]; then
        stop_reason="could not parse metrics at ${mb} Mb"; break; fi
    if [ "$ao" = "$na" ] && [ "$bo" = "$nb" ]; then ok="PASS"; else ok="FAIL"; fi

    now=$(date +%s); elapsed=$(( now - start ))
    printf '  %7s %11s %11s %9s %8s %7s %5s %10s %9s %9s\n' \
           "$mb" "$nA" "$nB" "$sp" "$ao" "$ok" "$dec" "$wall" "$mem" "$elapsed"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
           "$mb" "$bp" "$nA" "$nB" "$sp" "$ao" "$bo" "$ok" "$dec" "$wall" "$mem" >> "$TSV"

    if [ "$ok" = "FAIL" ]; then stop_reason="INCORRECT diff at ${mb} Mb - stopping"; break; fi

    ceiling="$mb Mb  |  $nA + $nB nodes  |  $sp bubbles  |  ${wall}s  |  ${mem} MB peak"
    last_nodes=$(( nA + nB )); last_mem="$mem"; last_wall="$wall"; last_mb="$mb"
    mb=$(( mb * 2 ))
done

echo
note "results table: $TSV"
note "stopped because: ${stop_reason:-done}"
if [ -n "$ceiling" ]; then
    echo "CEILING reached on this node:"
    echo "  $ceiling"
    awk -F'\t' 'NR>1 && $8=="PASS"{n=$3+$4; v=$2/'"$DENSITY"'; w=$10; mb=$1}
        END{ if(w>0) printf "  throughput at ceiling: %.0f variants/s, %.0f nodes/s (single core)\n", v/w, n/w }' "$TSV"
    echo "OK - pushed to the limit; diff stayed EXACT at every completed size"
else
    echo "No point completed - lower START_MB or check resources." >&2; exit 1
fi
