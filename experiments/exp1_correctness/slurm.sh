#!/bin/bash --login
# ============================================================
# SLURM job script — Experiment 1 (Correctness)
#
# Submit with:
#   sbatch experiments/exp1_correctness/slurm.sh
#
# ADJUST: --partition and --account for your HPC cluster.
# ============================================================
#SBATCH --job-name=g2g_exp1_correctness
#SBATCH --output=logs/exp1_%j.out
#SBATCH --error=logs/exp1_%j.err
#SBATCH --partition=mc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shakleyaahmed@gmail.com

set -euo pipefail

G2G_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "$G2G_ROOT/logs"

echo "[slurm] Job $SLURM_JOB_ID started on $(hostname) at $(date)"
echo "[slurm] CPUs: $SLURM_CPUS_ON_NODE  RAM: $(free -h | awk '/Mem/ {print $2}')"

# Load conda
source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || \
    source "$HOME/miniconda/etc/profile.d/conda.sh"
conda activate g2g_eval

cd "$G2G_ROOT"

# Record start time
echo "exp1_start_epoch=$(date +%s)" >> "results/phase5/exp1/slurm_timing_${SLURM_JOB_ID}.txt"

bash experiments/exp1_correctness/run.sh

echo "exp1_end_epoch=$(date +%s)" >> "results/phase5/exp1/slurm_timing_${SLURM_JOB_ID}.txt"
echo "[slurm] Exp1 completed at $(date)"
