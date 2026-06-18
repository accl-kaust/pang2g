#!/bin/bash --login
# SLURM — HPRC assembly download (Stage 1 of Exp6)
#SBATCH --job-name=g2g_exp6_download
#SBATCH --output=logs/exp6_download_%j.out
#SBATCH --error=logs/exp6_download_%j.err
#SBATCH --partition=mc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
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
export N_PER_POP=20

echo "[slurm] Download job $SLURM_JOB_ID started at $(date)"
bash experiments/download/download_hprc_assemblies.sh
echo "[slurm] Download completed at $(date)"
