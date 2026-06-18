#!/usr/bin/env bash
# SLURM — Experiment 2 (Minimap2 Consistency)
#SBATCH --job-name=g2g_exp2_minimap2
#SBATCH --output=logs/exp2_%j.out
#SBATCH --error=logs/exp2_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shakleyaahmed@gmail.com

set -euo pipefail
G2G_ROOT="$HOME/Research/G2G"
mkdir -p "$G2G_ROOT/logs"
source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || \
    source "$HOME/miniconda/etc/profile.d/conda.sh"
conda activate g2g_eval
cd "$G2G_ROOT"
echo "[slurm] Exp2 job $SLURM_JOB_ID started at $(date)"
bash experiments/exp2_minimap2/run.sh
echo "[slurm] Exp2 completed at $(date)"
