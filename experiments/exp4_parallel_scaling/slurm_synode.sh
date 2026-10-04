#!/bin/bash --login
# ============================================================
# SLURM — Experiment 4 (Parallel Scaling) — synode variant
#
# Targets the synode partition (384 CPUs / 1.5 TB RAM). This
# lets the thread sweep extend to 256 threads so the Amdahl
# fit is anchored well past the 128-core ceiling of the mc nodes.
#
# Submit: sbatch experiments/exp4_parallel_scaling/slurm_synode.sh
# ============================================================
#SBATCH --job-name=g2g_exp4_parallel
#SBATCH --output=logs/exp4_%j.out
#SBATCH --error=logs/exp4_%j.err
#SBATCH --partition=batch      # adjust: a 128+ core node (was a dedicated 384-core partition)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=768G
#SBATCH --time=14:00:00

set -euo pipefail

G2G_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "$G2G_ROOT/logs"

source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || \
    source "$HOME/miniconda/etc/profile.d/conda.sh"
conda activate g2g_eval

cd "$G2G_ROOT"

export G2G_DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
export MAX_THREADS="$SLURM_CPUS_ON_NODE"
export N_SNARLS="100000"
export REPEATS="5"

echo "[slurm] Job $SLURM_JOB_ID  CPUs=$SLURM_CPUS_ON_NODE  date=$(date)"

bash experiments/exp4_parallel_scaling/run.sh

echo "[slurm] Exp4 completed at $(date)"
