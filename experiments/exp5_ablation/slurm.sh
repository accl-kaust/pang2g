#!/usr/bin/env bash
# ============================================================
# SLURM — Experiment 5 (Ablation Study)
# Needs ~64 CPUs so parallel variants are stressed properly.
#
# Submit: sbatch experiments/exp5_ablation/slurm.sh
# ============================================================
#SBATCH --job-name=g2g_exp5_ablation
#SBATCH --output=logs/exp5_%j.out
#SBATCH --error=logs/exp5_%j.err
#SBATCH --partition=compute          # <-- change to your partition
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shakleyaahmed@gmail.com

set -euo pipefail

G2G_ROOT="$HOME/Research/G2G"
mkdir -p "$G2G_ROOT/logs"

source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || \
    source "$HOME/miniconda/etc/profile.d/conda.sh"
conda activate g2g_eval

cd "$G2G_ROOT"
export G2G_DATA_ROOT="/data/fat/g2g"
export THREADS="$SLURM_CPUS_ON_NODE"
export REPEATS="5"

echo "[slurm] Job $SLURM_JOB_ID  CPUs=$SLURM_CPUS_ON_NODE  date=$(date)"

bash experiments/exp5_ablation/run.sh

echo "[slurm] Exp5 completed at $(date)"
