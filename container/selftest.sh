#!/usr/bin/env bash
# Run the README's test case inside the container image.
#
#   docker run --rm g2g /opt/g2g/container/selftest.sh [--quick]
#   apptainer exec g2g.sif /opt/g2g/container/selftest.sh [--quick]
#   singularity exec g2g.sif /opt/g2g/container/selftest.sh [--quick]
#
# The image is read-only under Apptainer/Singularity, and the demos write
# their outputs next to the scripts, so the tree is copied to a writable
# scratch directory first and the pre-built binary in the image is used.
# Set G2G_SELFTEST_DIR to choose that directory (default: under $TMPDIR).
set -uo pipefail
SRC=/opt/g2g
WORK="${G2G_SELFTEST_DIR:-${TMPDIR:-/tmp}/g2g_selftest_$$}"
mkdir -p "$WORK" && cp -r "$SRC/." "$WORK/" && cd "$WORK" || exit 1
echo "selftest: image binary $SRC/build/release/g2g, work dir $WORK"
BUILD_DIR="$SRC/build/release" bash run_test_case.sh --no-build "$@"
