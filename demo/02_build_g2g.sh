#!/usr/bin/env bash
# Build the G2G binary (Release) at $REPO/build/g2g.
source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

log "Building G2G (Release)"
cmake -B "$REPO_DIR/build" -DCMAKE_BUILD_TYPE=Release -G Ninja "$REPO_DIR"
cmake --build "$REPO_DIR/build" -j"$(nproc)" --target g2g

"$G2G" --help 2>&1 | head -3 || true
echo "OK — built $G2G"
