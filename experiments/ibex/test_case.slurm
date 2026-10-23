#!/bin/bash --login
# PanG2G test case as a KAUST IBEX batch job.  From the repository root:
#
#     sbatch experiments/ibex/test_case.slurm
#     tail -f logs/test_case_<jobid>.out
#
# Loads the modules used for the paper's IBEX runs, builds the engine into
# build/ibex with -march=x86-64-v3 (so the same binary also runs on the login
# node and on every other IBEX node type), then runs run_test_case.sh: the
# planted ground-truth diff, the enrichment self-test and the synthetic Q4
# pipeline.  The GoogleTest suite is included when googletest is cached under
# /ibex/user (compute nodes have no outbound network; fetch it once, on the
# login node):
#
#     mkdir -p /ibex/user/$USER
#     curl -fsSL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz | tar xz -C /ibex/user/$USER
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
# The versions below are the ones the paper's IBEX jobs used (October 2026);
# `module avail gcc cmake python` lists what is current if one has been retired.
module load gcc/12.2.0 cmake/3.28.4/gnu-11.3.1 python/3.11.0 || module load gcc cmake python
cd "${SLURM_SUBMIT_DIR:?submit this script with sbatch from the repository root}"
echo "[ibex] $(date)  job=${SLURM_JOB_ID}  host=$(hostname)  $(lscpu | sed -n 's/^Model name: *//p')  cpus=$(nproc)"

UNIT=""
GT="/ibex/user/$USER/googletest-1.14.0"
if [ -d "$GT" ]; then
    export G2G_GTEST_DIR="$GT"; UNIT="--unit"
else
    echo "[ibex] no cached googletest at $GT - unit tests skipped (see the header of this file)"
fi

BUILD_DIR=build/ibex G2G_MARCH=x86-64-v3 JOBS="${SLURM_CPUS_PER_TASK:-8}" bash run_test_case.sh $UNIT
rc=$?
echo "[ibex] $(date)  exit=$rc"
exit $rc
