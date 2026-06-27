#!/usr/bin/env bash
# Demo 5 — ablation of G2G's design choices.
#
# Builds the fully-optimized baseline and one binary per ablation (each with a
# single optimization DISABLED via a CMake flag), then compares them at a
# CONTROL size (where the choice should not matter) and a TREATMENT size (where
# it should), reporting the slowdown (disabled / baseline) of the phase each
# choice affects. The expected pattern: control ~1.0x, treatment > 1.0x.
#
# Honest notes (verified against the code):
#   * The CMake option names are paper-facing; the build maps each to the
#     internal macro the optimization actually reads (see CMakeLists.txt).
#   * vEB layout (G2G_BFS_LAYOUT) is a no-op on the current depth-1 snarls, so
#     its rows should read ~1.0x — that is a real result, not a measurement bug.
#   * Static work-partition (G2G_STATIC_PARTITION) is not implemented in the
#     runtime, so it is omitted.
#
# Env: PTHREADS  threads for the scheduler ablations (default min(nproc,32))
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
ABL="$WORK/ablation"; mkdir -p "$ABL"
NPROC="$(nproc)"; PT="${PTHREADS:-$(( NPROC < 32 ? NPROC : 32 ))}"
if ninja --version >/dev/null 2>&1; then GEN="Ninja"; else GEN="Unix Makefiles"; fi

gen() { python3 "$DEMO_DIR/make_synthetic_pair.py" --chop-a 25 --chop-b 25 "$@" >/dev/null; }

build() {  # $1 = CMake flag ("" = baseline) -> echoes binary path
    local key="${1:-baseline}" d="$ABL/b_${1:-baseline}"
    local fl=(); [ -n "$1" ] && fl=("-D$1=ON")
    if [ ! -x "$d/g2g" ]; then
        note "  building $key …"
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

log "Ablation — each optimization OFF vs the optimized baseline"
note "building baseline + ablation binaries into $ABL (a few minutes)…"
BASE=$(build "")

note "generating control + treatment graphs…"
gen --outdir "$ABL/g_small"   --length 60000  --n-shared 400  --n-a 20 --n-b 20
gen --outdir "$ABL/g_bignode" --length 250000 --n-shared 400  --n-a 20 --n-b 20
gen --outdir "$ABL/g_k2"      --length 200000 --n-shared 3000 --n-a 30 --n-b 30
gen --outdir "$ABL/g_k7"      --length 200000 --n-shared 2000 --n-a 30 --n-b 30 --max-k 7  --complex-frac 1.0 --alt-len 20
gen --outdir "$ABL/g_sv_s"    --length 100000 --n-shared 400  --n-a 20 --n-b 20 --max-k 40 --complex-frac 1.0 --alt-len 200
gen --outdir "$ABL/g_sv_b"    --length 400000 --n-shared 3000 --n-a 50 --n-b 50 --max-k 40 --complex-frac 1.0 --alt-len 200

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
row "vEB layout"              G2G_BFS_LAYOUT     dp    "--threads 1"    "$ABL/g_sv_s"  "$ABL/g_sv_b"

echo
note "slowdown = disabled / baseline   (1.00x = no effect; >1 = the optimization helps here)"
note "results: $ABL/ablation_results.tsv"
note "vEB rows read ~1.0x: the layout is a no-op on depth-1 snarls (honest null)."
note "static work-partition is not implemented in the runtime — omitted."
echo "OK — ablation complete"
