#!/bin/bash --login
# ============================================================
# SLURM - Experiment 3 (Throughput Scaling) - large-node variant
#
# Same as slurm.sh but targets a node with more cores
# (384 CPUs / 1.5 TB RAM per node). Gives the largest graphs
# (500K–1M snarls) more headroom and faster parallel runs.
#
# Submit with:
#   sbatch experiments/exp3_throughput/slurm_large_node.sh
# ============================================================
#SBATCH --job-name=g2g_exp3_throughput
#SBATCH --output=logs/exp3_%j.out
#SBATCH --error=logs/exp3_%j.err
#SBATCH --partition=batch      # adjust: a 128+ core node (was a dedicated 384-core partition)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=512G
#SBATCH --time=14:00:00

set -euo pipefail

G2G_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "$G2G_ROOT/logs"

source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || \
    source "$HOME/miniconda/etc/profile.d/conda.sh"
conda activate g2g_eval

cd "$G2G_ROOT"

export G2G_DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
export THREADS="$SLURM_CPUS_ON_NODE"

echo "[slurm] Job $SLURM_JOB_ID  node=$(hostname)  CPUs=$SLURM_CPUS_ON_NODE  date=$(date)"
mkdir -p "results/phase5/exp3"
echo "exp3_start_epoch=$(date +%s)" \
    >> "results/phase5/exp3/slurm_timing_${SLURM_JOB_ID}.txt"

bash experiments/exp3_throughput/run.sh

echo "exp3_end_epoch=$(date +%s)" \
    >> "results/phase5/exp3/slurm_timing_${SLURM_JOB_ID}.txt"
echo "[slurm] Exp3 completed at $(date)"
