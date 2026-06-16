#!/bin/bash
# =============================================================================
# lib_measure.sh — Shared measurement & logging library for G2G Phase 5
# =============================================================================
# Source this from every experiment script:
#     source "$(dirname "$0")/lib_measure.sh"
#
# Provides:
#   measure_full      — run a command, capture wall / peak-RSS / CPU% / exit code
#   wall_only         — wall-clock seconds only (fast path)
#   l3_miss_rate      — perf-based L3 miss rate (gracefully degrades to N/A)
#   capture_env       — dump full hardware + software environment to a file
#   log               — timestamped stderr log line
#   tsv_init          — create a TSV with a header if it does not exist
#
# All functions are pure bash + coreutils + (optional) perf. No Python needed
# for the measurement path itself, so a broken conda env cannot corrupt timings.
#
# Design notes for reproducibility (IEEE BIBM):
#   * Wall clock is taken from GNU /usr/bin/time -v, NOT bash `time`, because
#     bash `time` cannot report peak RSS and rounds to 1ms.
#   * Peak RSS is "Maximum resident set size (kbytes)" — the high-water mark of
#     physical memory, which is the number reviewers expect for "memory used".
#   * Every measured command also records User+System CPU seconds so that
#     algorithmic cost (CPU-seconds) and wall-clock cost are never conflated —
#     this is what lets us reason about parallel efficiency honestly.
# =============================================================================

# ── Locate GNU time ──────────────────────────────────────────────────────────
# Some clusters ship it as /usr/bin/time, some as `gtime` (conda), some only as
# the bash builtin (useless to us). Resolve once.
_resolve_gnu_time() {
    if [[ -x /usr/bin/time ]]; then
        GNU_TIME=/usr/bin/time
    elif command -v gtime >/dev/null 2>&1; then
        GNU_TIME=$(command -v gtime)
    else
        GNU_TIME=""
    fi
    export GNU_TIME
}
_resolve_gnu_time

# ── Timestamped logging ──────────────────────────────────────────────────────
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

# ── TSV initialisation (idempotent) ──────────────────────────────────────────
# tsv_init <path> <tab-separated-header>
tsv_init() {
    local path="$1"; shift
    local header="$1"
    if [[ ! -f "$path" ]]; then
        printf '%b\n' "$header" > "$path"
    fi
}

# ── Full measurement ─────────────────────────────────────────────────────────
# measure_full <command...>
# Echoes a single line:  <wall_sec> <peak_rss_kb> <user_sec> <sys_sec> <cpu_pct> <exit_code>
# On failure (non-zero exit), still returns the captured metrics with the real
# exit code so the caller can decide whether to record or abort.
measure_full() {
    measure_full_capture /dev/null "$@"
}

# ── Full measurement, preserving the command's own stderr ────────────────────
# measure_full_capture <stderr_dest> <command...>
# Same as measure_full but the COMMAND's stderr is written to <stderr_dest> so
# the caller can scrape application log lines (e.g. G2G's "[timing] ..." line).
# GNU time's own report is kept separate in an internal tempfile.
#
# Implementation detail: GNU time writes its report to the fd given by -o. We
# point -o at a private tempfile and let the command's stderr flow to the
# caller's destination. This is the ONLY correct way to get both streams.
measure_full_capture() {
    local stderr_dest="$1"; shift
    local rep; rep=$(mktemp)
    local exit_code=0

    if [[ -n "$GNU_TIME" ]]; then
        # -o rep  : time report -> private tempfile
        # 2>dest  : command's stderr -> caller destination
        "$GNU_TIME" -v -o "$rep" "$@" 2>"$stderr_dest" || exit_code=$?
        local wall peak user sys cpu
        wall=$(grep "Elapsed (wall clock)" "$rep" | awk '{print $NF}' | \
               awk -F: '{if (NF==3) print $1*3600+$2*60+$3; else print $1*60+$2}')
        peak=$(grep "Maximum resident set size" "$rep" | awk '{print $NF}')
        user=$(grep "User time (seconds)"   "$rep" | awk '{print $NF}')
        sys=$( grep "System time (seconds)" "$rep" | awk '{print $NF}')
        cpu=$( grep "Percent of CPU this job got" "$rep" | awk '{print $NF}' | tr -d '%')
        wall=${wall:-0}; peak=${peak:-0}; user=${user:-0}; sys=${sys:-0}; cpu=${cpu:-0}
        echo "$wall $peak $user $sys $cpu $exit_code"
    else
        # Fallback: SECONDS for wall, no peak RSS. Preserve cmd stderr to dest.
        local start end wall
        start=$(date +%s.%N)
        "$@" >/dev/null 2>"$stderr_dest" || exit_code=$?
        end=$(date +%s.%N)
        wall=$(awk "BEGIN{print $end-$start}")
        echo "$wall 0 0 0 0 $exit_code"
    fi
    rm -f "$rep"
    return 0
}

# ── Wall-only (fast path for inner sweep loops) ──────────────────────────────
wall_only() {
    local tmp; tmp=$(mktemp)
    if [[ -n "$GNU_TIME" ]]; then
        "$GNU_TIME" -v "$@" 2>"$tmp" || true
        grep "Elapsed (wall clock)" "$tmp" | awk '{print $NF}' | \
            awk -F: '{if (NF==3) print $1*3600+$2*60+$3; else print $1*60+$2}'
    else
        local start end
        start=$(date +%s.%N); "$@" >/dev/null 2>&1 || true; end=$(date +%s.%N)
        awk "BEGIN{print $end-$start}"
    fi
    rm -f "$tmp"
}

# ── perf L3 miss rate (graceful) ─────────────────────────────────────────────
HAVE_PERF=0
if command -v perf >/dev/null 2>&1 && perf stat -e cache-references true >/dev/null 2>&1; then
    HAVE_PERF=1
fi
export HAVE_PERF

l3_miss_rate() {
    if [[ $HAVE_PERF -eq 0 ]]; then echo "NA"; return; fi
    local tmp; tmp=$(mktemp)
    perf stat -e cache-references,cache-misses "$@" 2>"$tmp" >/dev/null || true
    local refs misses
    refs=$(  grep "cache-references" "$tmp" | head -1 | awk '{gsub(",","",$1); print $1}')
    misses=$(grep "cache-misses"     "$tmp" | head -1 | awk '{gsub(",","",$1); print $1}')
    rm -f "$tmp"
    if [[ -z "$refs" || "$refs" == "0" || -z "$misses" ]]; then
        echo "NA"
    else
        awk "BEGIN{printf \"%.4f\", $misses/$refs}"
    fi
}

# ── Environment capture for reproducibility appendix ─────────────────────────
# capture_env <output_file>
capture_env() {
    local out="$1"
    {
        echo "# G2G Phase 5 — Run Environment"
        echo "# Generated: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        echo ""
        echo "## Host"
        echo "hostname: $(hostname)"
        echo "kernel:   $(uname -r)"
        echo "os:       $(. /etc/os-release 2>/dev/null && echo "$PRETTY_NAME" || echo unknown)"
        echo ""
        echo "## CPU"
        lscpu 2>/dev/null | grep -E "Model name|^CPU\(s\)|Thread|Core|Socket|NUMA node\(s\)|L1d|L1i|L2|L3|MHz" || true
        echo ""
        echo "## Memory"
        free -h 2>/dev/null | head -2 || true
        echo ""
        echo "## NUMA topology"
        numactl --hardware 2>/dev/null | head -20 || echo "numactl not available"
        echo ""
        echo "## Toolchain"
        echo "g++:      $(g++ --version 2>/dev/null | head -1 || echo missing)"
        echo "cmake:    $(cmake --version 2>/dev/null | head -1 || echo missing)"
        echo "openmp:   $(echo | g++ -fopenmp -dM -E - 2>/dev/null | grep -i _OPENMP || echo 'n/a')"
        echo "minimap2: $(minimap2 --version 2>/dev/null || echo missing)"
        echo "perf:     $(perf --version 2>/dev/null || echo 'missing / restricted')"
        echo ""
        echo "## SLURM"
        echo "job_id:        ${SLURM_JOB_ID:-NA}"
        echo "job_name:      ${SLURM_JOB_NAME:-NA}"
        echo "node:          ${SLURMD_NODENAME:-NA}"
        echo "cpus_per_task: ${SLURM_CPUS_PER_TASK:-NA}"
        echo "mem_per_node:  ${SLURM_MEM_PER_NODE:-NA}"
        echo "partition:     ${SLURM_JOB_PARTITION:-NA}"
        echo ""
        echo "## perf_event_paranoid"
        cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo unknown
    } > "$out"
    log "environment captured -> $out"
}

# ── SHA-256 of a file (reproducibility manifest) ─────────────────────────────
sha_of() {
    sha256sum "$1" 2>/dev/null | awk '{print $1}' || echo "NA"
}
