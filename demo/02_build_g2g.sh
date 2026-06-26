#!/usr/bin/env bash
# Build the G2G binary (Release) at $REPO/build/g2g.
#
# Robust to two common breakages:
#   * a build/ cache pinned to a build tool that no longer exists
#     (e.g. a removed conda 'ninja') — we detect and wipe it;
#   * no working 'ninja' on PATH — we fall back to Unix Makefiles.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

BUILD="$REPO_DIR/build"

# 1. If a previous configure recorded a make program that is now missing
#    (the classic ".../envs/<env>/bin/ninja: no such file" error), the cache is
#    unusable — wipe it so we can reconfigure from scratch.
if [ -f "$BUILD/CMakeCache.txt" ]; then
    mk=$(sed -n 's/^CMAKE_MAKE_PROGRAM:[^=]*=//p' "$BUILD/CMakeCache.txt" | head -1)
    if [ -n "$mk" ] && [ ! -x "$mk" ]; then
        note "build cache points at a missing build tool ($mk) — wiping build/"
        rm -rf "$BUILD"
    fi
fi

# 2. Choose a generator that actually works.
if ninja --version >/dev/null 2>&1; then
    GEN="Ninja"
else
    GEN="Unix Makefiles"
fi
log "Building G2G (Release, generator: $GEN)"

# -DG2G_BUILD_TESTS=OFF: the demo only needs the g2g binary, so skip the
# googletest fetch entirely (no network needed for the build).
configure() { cmake -B "$BUILD" -G "$GEN" -DCMAKE_BUILD_TYPE=Release \
                    -DG2G_BUILD_TESTS=OFF "$REPO_DIR"; }
# 3. If configuration still fails (e.g. a generator mismatch in an old cache),
#    wipe and try once more from a clean directory.
configure || { note "reconfiguring build/ from scratch"; rm -rf "$BUILD"; configure; }

cmake --build "$BUILD" -j"$(nproc)" --target g2g

"$G2G" --help 2>&1 | head -3 || true
echo "OK — built $G2G"
