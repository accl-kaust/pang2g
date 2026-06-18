#!/bin/bash --login
# ============================================================
# SLURM — Experiment 3 (Throughput Scaling) — synode variant
#
# Same as slurm.sh but targets the high-core synode partition
# (384 CPUs / 1.5 TB RAM per node). Gives the largest graphs
# (500K–1M snarls) more headroom and faster parallel runs.
#
# Submit with:
#   sbatch experiments/exp3_throughput/slurm_synode.sh
# ============================================================
#SBATCH --job-name=g2g_exp3_throughput
#SBATCH --output=logs/exp3_%j.out
#SBATCH --error=logs/exp3_%j.err
#SBATCH --partition=synode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=512G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shakleyaahmed@gmail.com

set -euo pipefail

G2G_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "$G2G_ROOT/logs"

source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || \
    source "$HOME/miniconda/etc/profile.d/conda.sh"
conda activate g2g_eval

cd "$G2G_ROOT"

export G2G_DATA_ROOT="/data/fat/g2g"
export THREADS="$SLURM_CPUS_ON_NODE"

echo "[slurm] Job $SLURM_JOB_ID  node=$(hostname)  CPUs=$SLURM_CPUS_ON_NODE  date=$(date)"
mkdir -p "results/phase5/exp3"
echo "exp3_start_epoch=$(date +%s)" \
    >> "results/phase5/exp3/slurm_timing_${SLURM_JOB_ID}.txt"

bash experiments/exp3_throughput/run.sh

echo "exp3_end_epoch=$(date +%s)" \
    >> "results/phase5/exp3/slurm_timing_${SLURM_JOB_ID}.txt"
echo "[slurm] Exp3 completed at $(date)"
