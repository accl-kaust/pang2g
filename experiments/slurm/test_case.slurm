#!/bin/bash --login
# PanG2G test case as a Slurm batch job.  From the repository root:
#
#     sbatch experiments/slurm/test_case.slurm
#     tail -f logs/test_case_<jobid>.out
#
# Adjust --partition (and add --account) for your site.  If your cluster
# provides its toolchain through environment modules, name them in
# G2G_MODULES; the job loads them before building:
#
#     G2G_MODULES="gcc cmake python" sbatch experiments/slurm/test_case.slurm
#
# The job builds the engine into build/slurm with -march=x86-64-v3 (so the same
# binary also runs on the login node and on any other node type) and runs
# run_test_case.sh: the planted ground-truth diff, the enrichment self-test and
# the synthetic biological validation.  The GoogleTest suite is included when
# G2G_GTEST_DIR points at an unpacked googletest-1.14.0 source tree; compute
# nodes often have no outbound network, so fetch it once on the login node:
#
#     curl -fsSL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz | tar xz -C /path/to/scratch
#     G2G_GTEST_DIR=/path/to/scratch/googletest-1.14.0 sbatch experiments/slurm/test_case.slurm
#
# About 3 minutes on 8 cores; the compile dominates.
#SBATCH --job-name=g2g_test_case
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/test_case_%j.out
set -uo pipefail
if [ -n "${G2G_MODULES:-}" ]; then
    module load $G2G_MODULES || { echo "[slurm] module load $G2G_MODULES failed"; exit 1; }
fi
cd "${SLURM_SUBMIT_DIR:?submit this script with sbatch from the repository root}"
echo "[slurm] $(date)  job=${SLURM_JOB_ID}  host=$(hostname)  $(lscpu | sed -n 's/^Model name: *//p')  cpus=$(nproc)"

UNIT=""
if [ -n "${G2G_GTEST_DIR:-}" ] && [ -d "$G2G_GTEST_DIR" ]; then
    UNIT="--unit"
else
    echo "[slurm] G2G_GTEST_DIR not set: unit tests skipped (see the header of this file)"
fi

BUILD_DIR=build/slurm G2G_MARCH=x86-64-v3 JOBS="${SLURM_CPUS_PER_TASK:-8}" bash run_test_case.sh $UNIT
rc=$?
echo "[slurm] $(date)  exit=$rc"
exit $rc
