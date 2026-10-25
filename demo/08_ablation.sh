#!/usr/bin/env bash
# Demo 5 - ablation of G2G's design choices.
#
# Builds the fully-optimized baseline and one binary per ablation (each with a
# single optimization DISABLED via a CMake flag), then compares them at a
# CONTROL size (where the choice should not matter) and a TREATMENT size (where
# it should), reporting the slowdown (disabled / baseline) of the phase each
# choice affects. The expected pattern: control ~1.0x, treatment > 1.0x.
#
# Notes:
#   * The CMake option names are user-facing; the build maps each to the
#     internal macro the optimization actually reads (see CMakeLists.txt).
#   * vEB layout (G2G_BFS_LAYOUT) is a no-op on the current depth-1 snarls, so
#     its rows should read ~1.0x - that is a real result, not a measurement bug.
#   * Static work-partition (G2G_STATIC_PARTITION) replaces work-stealing with a
#     fixed cyclic (ti % P) task assignment. Its control is a BALANCED workload
#     (every complex snarl identical, via --uniform-k) where a static split is
#     already even (~1.0x); its treatment is a heavy-tailed workload where a few
#     threads draw the big snarls and cannot shed them, so static loses (>1.0x).
#
# Env: PTHREADS      threads for the steal-largest ablation (default min(nproc,32))
#      SCHED_THREADS threads for the static-partition ablation (default min(nproc,8)):
#                    work-stealing's edge over a static split is sharpest when the
#                    heavy-snarl count is a small multiple of P; far past that both
#                    leave threads idle on a finite heavy tail, so it is capped for
#                    a machine-independent result (raise it on a big node to taste).
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
ABL="$WORK/ablation"; mkdir -p "$ABL"
NPROC="$(nproc)"; PT="${PTHREADS:-$(( NPROC < 32 ? NPROC : 32 ))}"
ST="${SCHED_THREADS:-$(( NPROC < 8 ? NPROC : 8 ))}"
if ninja --version >/dev/null 2>&1; then GEN="Ninja"; else GEN="Unix Makefiles"; fi

gen() { python3 "$DEMO_DIR/make_synthetic_pair.py" --chop-a 25 --chop-b 25 "$@" >/dev/null; }

build() {  # $1 = CMake flag ("" = baseline) -> echoes binary path
    local key="${1:-baseline}" d="$ABL/b_${1:-baseline}"
    local fl=(); [ -n "$1" ] && fl=("-D$1=ON")
    if [ ! -x "$d/g2g" ]; then
        note "  building $key …" >&2    # progress to stderr; only the path on stdout
        cmake -B "$d" -G "$GEN" -DCMAKE_BUILD_TYPE=Release -DG2G_BUILD_TESTS=OFF \
              "${fl[@]}" "$REPO_DIR" >/dev/null 2>&1
        cmake --build "$d" -j"$NPROC" --target g2g >/dev/null 2>&1
    fi
    echo "$d/g2g"
}

metric() {  # $1=bin $2=A $3=B $4=timing-field $5..=run flags -> min over 3 reps
    local bin="$1" A="$2" B="$3" f="$4"; shift 4
    local best="" v
    for _ in 1 2 3; do
        v=$(timeout 900 "$bin" "$A" "$B" "$@" --verbose-timing 2>&1 >/dev/null \
            | sed -n "s/.* ${f}=\([0-9.]*\).*/\1/p")
        [ -z "$v" ] && continue
        { [ -z "$best" ] || awk -v x="$v" -v y="$best" 'BEGIN{exit !(x<y)}'; } && best="$v"
    done
    echo "${best:-NA}"
}
ratio() { awk -v a="$1" -v b="$2" 'BEGIN{ if (b>0 && a!="NA") printf "%.2fx", a/b; else printf "NA" }'; }

log "Ablation - each optimization OFF vs the optimized baseline"
note "building baseline + ablation binaries into $ABL (a few minutes)…"
BASE=$(build "")

note "generating control + treatment graphs…"
# Anchor control is sized BELOW the minimiser's internal threshold (≈1e6 naïve
# comparisons, i.e. ~1000 backbone nodes) so the baseline also uses the naïve
# scan there → a true ~1.0× control; the treatment is well above it.
gen --outdir "$ABL/g_small"   --length 18000  --n-shared 200  --n-a 20 --n-b 20
gen --outdir "$ABL/g_bignode" --length 150000 --n-shared 200  --n-a 20 --n-b 20
gen --outdir "$ABL/g_k2"      --length 200000 --n-shared 3000 --n-a 30 --n-b 30
gen --outdir "$ABL/g_k7"      --length 200000 --n-shared 2000 --n-a 30 --n-b 30 --max-k 7  --complex-frac 1.0 --alt-len 20
gen --outdir "$ABL/g_sv_s"    --length 50000  --n-shared 120  --n-a 20 --n-b 20 --max-k 40 --complex-frac 1.0 --alt-len 200
gen --outdir "$ABL/g_sv_b"    --length 400000 --n-shared 3000 --n-a 50 --n-b 50 --max-k 40 --complex-frac 1.0 --alt-len 200
# static-partition: BALANCED control (identical heavy snarls) vs SKEWED treatment
# (heavy-tailed) - same idea, differing only in per-snarl work variance.
gen --outdir "$ABL/g_bal"     --length 250000 --n-shared 200  --n-a 10 --n-b 10 --max-k 16 --complex-frac 1.0 --alt-len 120 --uniform-k
gen --outdir "$ABL/g_skew"    --length 400000 --n-shared 320  --n-a 10 --n-b 10 --max-k 40 --complex-frac 1.0 --alt-len 160

row() {  # $1=label $2=flag $3=metric $4=runflags $5=ctl-dir $6=trt-dir
    local bin; bin=$(build "$2")
    local cb cd tb td
    cb=$(metric "$BASE" "$5/demoA.gfa" "$5/demoB.gfa" "$3" $4)
    cd=$(metric "$bin"  "$5/demoA.gfa" "$5/demoB.gfa" "$3" $4)
    tb=$(metric "$BASE" "$6/demoA.gfa" "$6/demoB.gfa" "$3" $4)
    td=$(metric "$bin"  "$6/demoA.gfa" "$6/demoB.gfa" "$3" $4)
    printf '  %-26s %-9s %10s %12s\n' "$1" "$3" "$(ratio "$cd" "$cb")" "$(ratio "$td" "$tb")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$cb" "$cd" "$tb" "$td" >> "$ABL/ablation_results.tsv"
}

printf 'label\tflag\tphase\tctrl_base\tctrl_disabled\ttrt_base\ttrt_disabled\n' > "$ABL/ablation_results.tsv"
echo
printf '  %-26s %-9s %10s %12s\n' "optimization (disabled)" "phase" "control" "treatment"
printf '  %-26s %-9s %10s %12s\n' "-----------------------" "-----" "-------" "---------"
row "LAPJV matching"          G2G_BRUTE_MATCHING dp    "--threads 1"    "$ABL/g_k2"    "$ABL/g_k7"
row "minimiser anchor index"  G2G_NAIVE_ANCHOR   index "--threads 1"    "$ABL/g_small" "$ABL/g_bignode"
row "steal-largest policy"    G2G_STEAL_RANDOM   dp    "--threads $PT"  "$ABL/g_sv_s"  "$ABL/g_sv_b"
row "static work partition"   G2G_STATIC_PARTITION dp  "--threads $ST"  "$ABL/g_bal"   "$ABL/g_skew"
row "vEB layout"              G2G_BFS_LAYOUT     dp    "--threads 1"    "$ABL/g_sv_s"  "$ABL/g_sv_b"

echo
note "slowdown = disabled / baseline   (1.00x = no effect; >1 = the optimization helps here)"
note "results: $ABL/ablation_results.tsv"
note "vEB rows read ~1.0x: the layout is a no-op on depth-1 snarls, so 1.0x is the expected result."
note "static-partition control ~1.0x (balanced work) vs treatment >1.0x (skewed work)."
echo "OK - ablation complete"
