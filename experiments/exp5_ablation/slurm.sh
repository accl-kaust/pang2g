#!/bin/bash --login
# ============================================================
# SLURM - Experiment 5 (Ablation Study)
# Needs ~64 CPUs so parallel variants are stressed properly.
#
# Submit: sbatch experiments/exp5_ablation/slurm.sh
# ============================================================
#SBATCH --job-name=g2g_exp5_ablation
#SBATCH --output=logs/exp5_%j.out
#SBATCH --error=logs/exp5_%j.err
#SBATCH --partition=batch      # adjust to your cluster
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=12:00:00

set -euo pipefail

G2G_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "$G2G_ROOT/logs"

source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || \
    source "$HOME/miniconda/etc/profile.d/conda.sh"
conda activate g2g_eval

cd "$G2G_ROOT"
export G2G_DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
export THREADS="$SLURM_CPUS_ON_NODE"
export REPEATS="5"

echo "[slurm] Job $SLURM_JOB_ID  CPUs=$SLURM_CPUS_ON_NODE  date=$(date)"

bash experiments/exp5_ablation/run.sh

echo "[slurm] Exp5 completed at $(date)"
