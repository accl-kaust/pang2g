#!/usr/bin/env bash
# SLURM — Build HPRC pangenomes (Stage 2 of Exp6)
# ~48 CPU-hours per pangenome, 2 pangenomes = ~4–5 days
#SBATCH --job-name=g2g_exp6_build
#SBATCH --output=logs/exp6_build_%j.out
#SBATCH --error=logs/exp6_build_%j.err
#SBATCH --partition=compute          # <-- change to your partition
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=120:00:00
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

echo "[slurm] Build job $SLURM_JOB_ID started at $(date)"
bash experiments/exp6_hprc/02_build_pangenomes.sh
echo "[slurm] Build completed at $(date)"
