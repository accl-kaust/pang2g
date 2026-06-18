#!/usr/bin/env bash
# ============================================================
# Phase 5 — Build All G2G Variants
#
# Builds the production binary and five ablation variants.
# Each variant disables one optimization so we can measure
# its contribution in Experiment 5.
#
# Usage:
#   conda activate g2g_eval
#   bash experiments/01_build_all.sh
#
# Outputs (all under build/):
#   build/release/g2g                    production binary
#   build/abl_bfs_layout/g2g             ablation: BFS layout (no vEB)
#   build/abl_steal_random/g2g           ablation: random steal (no largest-first)
#   build/abl_naive_anchor/g2g           ablation: naive anchor (no minimizer index)
#   build/abl_brute_matching/g2g         ablation: brute-force matching (no Hungarian)
#   build/abl_static_partition/g2g       ablation: static partition (no work-stealing)
# ============================================================
set -euo pipefail

G2G_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NCPU=$(nproc)
ENV_NAME="${G2G_ENV_NAME:-g2g_eval}"

# ── Ensure the g2g_eval toolchain env is active ───────────────
# The cmake/ninja/compilers live inside the conda env, not in base.
if [[ "${CONDA_DEFAULT_ENV:-}" != "$ENV_NAME" ]]; then
    echo "[build] Activating conda env '$ENV_NAME' (was: ${CONDA_DEFAULT_ENV:-none})"
    # shellcheck disable=SC1090
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$ENV_NAME"
fi

# ── Preflight: confirm the build toolchain is on PATH ─────────
missing=()
for tool in cmake ninja; do
    command -v "$tool" &>/dev/null || missing+=("$tool")
done
# A C++ compiler from the env (cxx-compiler installs x86_64-conda-linux-gnu-g++)
if ! command -v c++ &>/dev/null && ! command -v g++ &>/dev/null \
   && ! ls "${CONDA_PREFIX:-/nonexistent}"/bin/*-g++ &>/dev/null; then
    missing+=("c++ compiler")
fi
if (( ${#missing[@]} > 0 )); then
    echo "[build] ERROR: missing build tools: ${missing[*]}" >&2
    echo "[build] Run 'bash experiments/00_setup_env.sh' first, then re-run this script." >&2
    exit 1
fi
echo "[build] Toolchain OK: cmake=$(command -v cmake)  ninja=$(command -v ninja)"

build_variant() {
    local name="$1"; shift
    local extra_flags=("$@")
    local dir="$G2G_ROOT/build/$name"
    echo ""
    echo "══════════════════════════════════════════"
    echo " Building: $name"
    echo "══════════════════════════════════════════"
    cmake -B "$dir" \
        -DCMAKE_BUILD_TYPE=Release \
        -G Ninja \
        "${extra_flags[@]}" \
        "$G2G_ROOT"
    cmake --build "$dir" -j"$NCPU"
    echo "[build] $name → $dir/g2g"
}

# Production build (all optimizations on)
build_variant release

# Ablation: disable cache-oblivious vEB layout → use BFS order
build_variant abl_bfs_layout        -DG2G_BFS_LAYOUT=ON

# Ablation: steal random task instead of largest
build_variant abl_steal_random      -DG2G_STEAL_RANDOM=ON

# Ablation: naive O(|VA||VB|) anchor scan instead of minimizer index
build_variant abl_naive_anchor      -DG2G_NAIVE_ANCHOR=ON

# Ablation: brute-force bipartite matching (exponential in k)
build_variant abl_brute_matching    -DG2G_BRUTE_MATCHING=ON

# Ablation: static round-robin partition instead of work-stealing
build_variant abl_static_partition  -DG2G_STATIC_PARTITION=ON

echo ""
echo "All builds complete. Summary:"
for variant in release abl_bfs_layout abl_steal_random abl_naive_anchor \
               abl_brute_matching abl_static_partition; do
    bin="$G2G_ROOT/build/$variant/g2g"
    if [[ -f "$bin" ]]; then
        size=$(du -sh "$bin" | cut -f1)
        echo "  ✓  $variant  ($size)"
    else
        echo "  ✗  $variant  MISSING"
    fi
done
