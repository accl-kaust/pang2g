#!/usr/bin/env bash
# ============================================================
# SLURM — Experiment 3 (Throughput Scaling)
#
# Runs all graph sizes (1K → 1M snarls) × 3 seeds × 4 tools.
# Estimated wall time: 8–12h depending on largest graphs.
#
# Submit with:
#   sbatch experiments/exp3_throughput/slurm.sh
#
# ADJUST: --partition, --account for your cluster.
# ============================================================
#SBATCH --job-name=g2g_exp3_throughput
#SBATCH --output=logs/exp3_%j.out
#SBATCH --error=logs/exp3_%j.err
#SBATCH --partition=compute          # <-- change to your partition
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=24:00:00
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

echo "[slurm] Job $SLURM_JOB_ID  node=$(hostname)  CPUs=$SLURM_CPUS_ON_NODE  date=$(date)"
echo "exp3_start_epoch=$(date +%s)" \
    >> "results/phase5/exp3/slurm_timing_${SLURM_JOB_ID}.txt"

bash experiments/exp3_throughput/run.sh

echo "exp3_end_epoch=$(date +%s)" \
    >> "results/phase5/exp3/slurm_timing_${SLURM_JOB_ID}.txt"
echo "[slurm] Exp3 completed at $(date)"
