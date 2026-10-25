#!/usr/bin/env bash
# run_test_case.sh - build PanG2G and run its end-to-end test case.
#
#   bash run_test_case.sh            # build + the three checks below        (~2 min)
#   bash run_test_case.sh --unit     # ... plus the GoogleTest suite (161 tests)
#   bash run_test_case.sh --real     # ... plus two real HLA-DRB1 pangenomes (downloads)
#
# Checks (each has a known expected answer; none needs network access or any
# third-party bioinformatics tool - only a C++ compiler, CMake and python3):
#
#   1. planted ground truth      demo/03_synthetic_demo.sh
#      Two graphs embed the same reference sequence but split it into different
#      nodes (25 bp vs 17 bp); 5 shared, 3 A-only and 2 B-only SNPs are planted.
#      `g2g --ref` must report exactly the 3 A-only and 2 B-only sites at the
#      planted coordinates, and a graph diffed against itself must report 0.
#   2. enrichment self-test      demo/hprc_enrich.py --self-test
#      The permutation statistic behind the biological validation must
#      recover a planted enrichment (fold >> 1, p < 0.01) and report none for a
#      random control (fold ~ 1).
#   3. biological validation, synthetic   demo/09_hprc_biology.sh --synthetic
#      A 10 Mb pair with 300 planted population-specific sites is diffed with
#      the real binary; the deltas must be enriched for a synthetic high-F_ST
#      gnomAD-SV fixture and for ClinVar-pathogenic sites, and a shuffled
#      control must not be.
#
# Options
#   --unit       configure with the GoogleTest suite and run ctest (CMake
#                downloads googletest; offline, set G2G_GTEST_DIR)
#   --real       download vg + minimap2 into demo/tools/ and diff two graphs
#                built from real HLA-DRB1 haplotypes (demo/04_real_hla_demo.sh)
#   --quick      skip check 3
#   --no-build   use the binary already in $BUILD_DIR
#
# Environment
#   BUILD_DIR      CMake build directory               (default: build/release)
#   G2G_MARCH      -march for the build: native (default), or x86-64-v3 to build
#                  on one node and run on another (login node -> compute node)
#   G2G_GTEST_DIR  unpacked googletest-1.14.0 source tree for an offline --unit
#   JOBS           parallel compile jobs               (default: nproc)
#
# Each step's output is also saved under logs/test_case/.  Exit status is 0
# only when every check passed.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

UNIT=0; REAL=0; QUICK=0; BUILD=1
for a in "$@"; do
    case "$a" in
        --unit) UNIT=1 ;;  --real) REAL=1 ;;  --quick) QUICK=1 ;;  --no-build) BUILD=0 ;;
        -h|--help) sed -n '2,/^set -uo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown option: $a (see --help)" >&2; exit 2 ;;
    esac
done
BUILD_DIR="${BUILD_DIR:-build/release}"
case "$BUILD_DIR" in /*) BUILD_ABS="$BUILD_DIR" ;; *) BUILD_ABS="$ROOT/$BUILD_DIR" ;; esac
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
LOGDIR="$ROOT/logs/test_case"; mkdir -p "$LOGDIR"

# -- toolchain ----------------------------------------------------------------
missing=()
command -v python3 >/dev/null 2>&1 || missing+=("python3")
CXX_BIN="${CXX:-}"
if [ "$BUILD" -eq 1 ]; then      # a compiler and CMake are only needed to build
    command -v cmake >/dev/null 2>&1 || missing+=("cmake")
    [ -n "$CXX_BIN" ] || CXX_BIN="$(command -v g++ || command -v c++ || command -v clang++ || true)"
    [ -n "$CXX_BIN" ] || missing+=("a C++ compiler (g++ / clang++)")
fi
if [ "${#missing[@]}" -gt 0 ]; then
    echo "ERROR: not found: ${missing[*]}" >&2
    echo "       on a cluster, load the toolchain modules first (e.g. module load gcc cmake python)" >&2
    exit 2
fi
echo "PanG2G test case  -  $(date)"
echo "  host: $(hostname)   cpus: $(nproc 2>/dev/null || echo '?')   python3: $(python3 --version 2>&1 | sed 's/Python //')"
if [ "$BUILD" -eq 1 ]; then
    echo "  compiler: $("$CXX_BIN" --version | head -1)   cmake: $(cmake --version | head -1 | sed 's/cmake version //')"
    echo "  build dir: $BUILD_DIR   -march: ${G2G_MARCH:-native}   jobs: $JOBS"
else
    echo "  binary: $BUILD_ABS/g2g (no build)"
fi

# -- step runner --------------------------------------------------------------
NAMES=(); STATUS=(); SECS=()
run_step() {   # <name> <command...>  - run, tee to a log, record PASS/FAIL
    local name="$1"; shift
    local slug t0 rc
    slug=$(printf '%s' "$name" | tr -c 'A-Za-z0-9' '_' | tr -s '_' | sed 's/_$//')
    t0=$(date +%s)
    printf '\n==== %s ====\n' "$name"
    "$@" 2>&1 | tee "$LOGDIR/$slug.log"
    rc=${PIPESTATUS[0]}
    NAMES+=("$name"); SECS+=("$(( $(date +%s) - t0 ))")
    if [ "$rc" -eq 0 ]; then STATUS+=("PASS")
    else STATUS+=("FAIL"); echo "  -> FAILED (exit $rc; log: $LOGDIR/$slug.log)"; fi
}
summary() {
    printf '\n==== summary ====\n'
    local i fail=0
    if [ "${#NAMES[@]}" -gt 0 ]; then
        for i in $(seq 0 $(( ${#NAMES[@]} - 1 ))); do
            printf '  %-4s  %-58s %5ss\n' "${STATUS[$i]}" "${NAMES[$i]}" "${SECS[$i]}"
            [ "${STATUS[$i]}" = FAIL ] && fail=1
        done
    fi
    echo
    if [ "${#NAMES[@]}" -eq 0 ]; then echo "NO CHECKS RAN"; return 1; fi
    if [ "$fail" -eq 0 ]; then echo "ALL CHECKS PASSED"; else echo "SOME CHECKS FAILED  (logs: $LOGDIR)"; fi
    return "$fail"
}

# -- build --------------------------------------------------------------------
build() {
    local gen="Unix Makefiles"; ninja --version >/dev/null 2>&1 && gen="Ninja"
    local tests=OFF; [ "$UNIT" -eq 1 ] && tests=ON
    local args=(-DCMAKE_BUILD_TYPE=Release "-DG2G_BUILD_TESTS=$tests")
    [ -n "${G2G_MARCH:-}" ]     && args+=("-DG2G_MARCH=$G2G_MARCH")
    [ -n "${G2G_GTEST_DIR:-}" ] && args+=("-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=$G2G_GTEST_DIR")
    # A cache whose build tool has vanished (e.g. a removed conda ninja) is unusable.
    if [ -f "$BUILD_ABS/CMakeCache.txt" ]; then
        local mk; mk=$(sed -n 's/^CMAKE_MAKE_PROGRAM:[^=]*=//p' "$BUILD_ABS/CMakeCache.txt" | head -1)
        if [ -n "$mk" ] && [ ! -x "$mk" ]; then echo "stale build cache ($mk missing): wiping $BUILD_DIR"; rm -rf "$BUILD_ABS"; fi
    fi
    echo "configure: cmake -B $BUILD_DIR -G \"$gen\" ${args[*]}"
    cmake -B "$BUILD_ABS" -G "$gen" "${args[@]}" "$ROOT" \
        || { echo "reconfiguring $BUILD_DIR from scratch"; rm -rf "$BUILD_ABS"; cmake -B "$BUILD_ABS" -G "$gen" "${args[@]}" "$ROOT"; } \
        || return 1
    if [ "$UNIT" -eq 1 ]; then cmake --build "$BUILD_ABS" -j"$JOBS"
    else                       cmake --build "$BUILD_ABS" -j"$JOBS" --target g2g; fi
}

export G2G="$BUILD_ABS/g2g"           # the demo scripts (demo/lib.sh) use this binary
[ "$BUILD" -eq 1 ] && run_step "build g2g ($BUILD_DIR)" build
if [ ! -x "$G2G" ]; then
    echo "ERROR: no binary at $G2G" >&2; summary; exit 1
fi
[ "$UNIT" -eq 1 ] && run_step "unit tests (ctest, $BUILD_DIR)" bash -c "cd '$BUILD_ABS' && ctest --output-on-failure -j'$JOBS'"

run_step "1. planted ground truth (demo/03_synthetic_demo.sh)"           bash demo/03_synthetic_demo.sh
run_step "2. enrichment self-test (demo/hprc_enrich.py --self-test)"     python3 demo/hprc_enrich.py --self-test
[ "$QUICK" -eq 0 ] && \
run_step "3. biological validation on synthetic data (demo/09_hprc_biology.sh --synthetic)" bash demo/09_hprc_biology.sh --synthetic
if [ "$REAL" -eq 1 ]; then
    run_step "real data: fetch vg + minimap2 (demo/01_setup_tools.sh)"    bash demo/01_setup_tools.sh
    run_step "real data: HLA-DRB1 pangenome diff (demo/04_real_hla_demo.sh)" bash demo/04_real_hla_demo.sh
fi
summary
