#!/usr/bin/env bash
# Shared configuration and helpers for the G2G demo scripts.
set -euo pipefail

DEMO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$DEMO_DIR/.." && pwd)"
TOOLS="$DEMO_DIR/tools"     # downloaded third-party binaries (git-ignored)
WORK="$DEMO_DIR/work"       # all generated data and outputs (git-ignored)
G2G="$REPO_DIR/build/g2g"

mkdir -p "$TOOLS" "$WORK"
export PATH="$TOOLS:$PATH"

if [ -t 1 ]; then
    log() { printf '\n\033[1m=== %s ===\033[0m\n' "$*"; }
else
    log() { printf '\n=== %s ===\n' "$*"; }   # plain when redirected to a log
fi
note() { printf '    %s\n' "$*"; }

# Resolve a tool by name: prefer demo/tools, then PATH. Empty if not found.
have() { command -v "$1" >/dev/null 2>&1; }
need_g2g() {
    [ -x "$G2G" ] || { echo "ERROR: $G2G not found — run demo/02_build_g2g.sh first" >&2; exit 1; }
}
