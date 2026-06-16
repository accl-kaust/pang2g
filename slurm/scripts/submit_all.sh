#!/bin/bash
# =============================================================================
# submit_all.sh — Master submission for the G2G Phase 5 evaluation campaign
# =============================================================================
# Submits the four experiment jobs to SLURM with the right ordering and shared
# graph cache, so the whole campaign can be launched with one command:
#
#     bash slurm/scripts/submit_all.sh
#
# ORDERING RATIONALE
#   * exp6 (biological) is the LONG POLE — HPRC download (~1 TB) + two
#     minigraph-cactus builds (~96 CPU-h). It shares nothing with the synthetic
#     experiments, so it is submitted FIRST and runs on its own node in
#     parallel with everything else.
#   * exp3 (throughput) GENERATES the synthetic graph cache that exp4 and exp5
#     reuse (via $G2G_GRAPHDIR). exp4 and exp5 therefore DEPEND on exp3 having
#     produced the graphs. We express this with `--dependency=afterok:<exp3>`.
#   * exp4 and exp5 are independent of each other and run concurrently once
#     exp3's graphs exist.
#
#   exp6 ─────────────────────────────────────────────►  (independent, long)
#   exp3 ──┬─► exp4                                       (afterok)
#          └─► exp5                                       (afterok)
#
# USAGE
#   bash submit_all.sh                 # submit everything
#   bash submit_all.sh --dry-run       # print sbatch commands, submit nothing
#   bash submit_all.sh exp3 exp4       # submit only the named jobs
#   bash submit_all.sh --no-bio        # skip exp6 (e.g. data not staged yet)
#   bash submit_all.sh --local         # run sequentially without SLURM (smoke)
#
# ENV OVERRIDES (passed through to the jobs via --export)
#   G2G_REPO, G2G_ENV, CONDA_BASE, G2G_GRAPHDIR, G2G_DATA_ROOT, G2G_SIZES,
#   G2G_RUNS, G2G_THREADS, G2G_N_TREATMENT, ...  (see each experiment script)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${G2G_REPO:-$(cd "$SCRIPT_DIR/../.." && pwd)}"
LOGDIR="$REPO/slurm/logs"
mkdir -p "$LOGDIR"

# Always export G2G_REPO so the submitted jobs can resolve repo paths even
# though SLURM runs the batch script from its spool dir (where BASH_SOURCE no
# longer points into the repo). --export=ALL below forwards this to every job.
export G2G_REPO="$REPO"

DRY_RUN=0
LOCAL=0
RUN_BIO=1
declare -a WHICH=()

for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        --local)   LOCAL=1 ;;
        --no-bio)  RUN_BIO=0 ;;
        exp3|exp4|exp5|exp6) WHICH+=("$arg") ;;
        -h|--help)
            grep '^#' "$0" | sed 's/^# \{0,1\}//' | head -40
            exit 0 ;;
        *) echo "unknown arg: $arg (use --help)"; exit 2 ;;
    esac
done

# Default selection = all four (minus bio if --no-bio).
if [[ ${#WHICH[@]} -eq 0 ]]; then
    WHICH=(exp3 exp4 exp5)
    [[ $RUN_BIO -eq 1 ]] && WHICH+=(exp6)
fi
want() { local x; for x in "${WHICH[@]}"; do [[ "$x" == "$1" ]] && return 0; done; return 1; }

EXP3="$SCRIPT_DIR/exp3_throughput.sh"
EXP4="$SCRIPT_DIR/exp4_parallel_scaling.sh"
EXP5="$SCRIPT_DIR/exp5_ablation.sh"
EXP6="$SCRIPT_DIR/exp6_biological.sh"

# Variables we forward into every job's environment.
EXPORT_VARS="ALL"
[[ -n "${G2G_EXPORT:-}" ]] && EXPORT_VARS="$G2G_EXPORT"

# ── Local mode: just run them in dependency order, no SLURM ──────────────────
if [[ $LOCAL -eq 1 ]]; then
    echo "[submit_all] LOCAL mode — running sequentially (no SLURM)"
    want exp6 && { echo "=== exp6 (bio) ==="; bash "$EXP6" || echo "exp6 returned $?"; }
    want exp3 && { echo "=== exp3 (throughput) ==="; bash "$EXP3"; }
    want exp4 && { echo "=== exp4 (parallel) ==="; bash "$EXP4"; }
    want exp5 && { echo "=== exp5 (ablation) ==="; bash "$EXP5"; }
    echo "[submit_all] LOCAL run complete."
    exit 0
fi

# ── SLURM mode ───────────────────────────────────────────────────────────────
# In dry-run we only print commands, so sbatch need not be installed.
if [[ $DRY_RUN -eq 0 ]]; then
    command -v sbatch >/dev/null 2>&1 || { echo "ERROR: sbatch not found and --local not set"; exit 1; }
fi

submit() {
    # submit <script> [extra sbatch args...] ; echoes the job id
    local script="$1"; shift
    local cmd=(sbatch --parsable --export="$EXPORT_VARS" "$@" "$script")
    if [[ $DRY_RUN -eq 1 ]]; then
        echo "[dry-run] ${cmd[*]}" >&2
        echo "DRYRUN"          # placeholder id
    else
        local id; id=$("${cmd[@]}")
        echo "$id"
    fi
}

declare -A JOBID

# exp6 first — independent, long-running.
if want exp6; then
    JOBID[exp6]=$(submit "$EXP6")
    echo "[submit_all] exp6 (biological)      -> ${JOBID[exp6]}"
fi

# exp3 next — produces the synthetic graph cache.
if want exp3; then
    JOBID[exp3]=$(submit "$EXP3")
    echo "[submit_all] exp3 (throughput)      -> ${JOBID[exp3]}"
fi

# exp4 / exp5 depend on exp3 (graphs). If exp3 wasn't submitted in this call,
# they run immediately (assuming the cache already exists on disk).
DEP=()
if want exp3 && [[ "${JOBID[exp3]}" != "DRYRUN" ]]; then
    DEP=(--dependency=afterok:"${JOBID[exp3]}")
fi

if want exp4; then
    JOBID[exp4]=$(submit "$EXP4" "${DEP[@]}")
    echo "[submit_all] exp4 (parallel scaling) -> ${JOBID[exp4]}  ${DEP[*]:-(no dep)}"
fi
if want exp5; then
    JOBID[exp5]=$(submit "$EXP5" "${DEP[@]}")
    echo "[submit_all] exp5 (ablation)         -> ${JOBID[exp5]}  ${DEP[*]:-(no dep)}"
fi

echo ""
if [[ $DRY_RUN -eq 1 ]]; then
    echo "[submit_all] dry-run complete — nothing submitted."
else
    echo "[submit_all] all jobs submitted. Monitor with:  squeue -u \$USER"
    echo "[submit_all] results land under: $REPO/results/{exp3,exp4,exp5,exp6}/"
fi
