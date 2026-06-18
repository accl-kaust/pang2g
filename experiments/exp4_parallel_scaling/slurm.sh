#!/bin/bash --login
# ============================================================
# SLURM — Experiment 4 (Parallel Scaling)
# Needs 128 CPUs and 256 GB RAM for HPRC-scale graph.
#
# Submit: sbatch experiments/exp4_parallel_scaling/slurm.sh
# ADJUST: --partition for your cluster.
# ============================================================
#SBATCH --job-name=g2g_exp4_parallel
#SBATCH --output=logs/exp4_%j.out
#SBATCH --error=logs/exp4_%j.err
#SBATCH --partition=mc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shakleyaahmed@gmail.com

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
