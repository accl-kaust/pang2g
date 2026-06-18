#!/bin/bash --login
# ============================================================
# DEPRECATED — do NOT submit this as a Slurm job.
#
# The compute nodes on this cluster have NO internet access, so the Exp6
# downloads cannot run in a batch job (they hang / fail immediately). Run the
# downloads on the LOGIN NODE inside tmux instead:
#
#   tmux new -s g2g_dl
#   conda activate g2g_eval
#   bash experiments/download/download_all_login.sh
#
# This wrapper is kept only to redirect anyone who still `sbatch`es it.
# ============================================================
#SBATCH --job-name=g2g_exp6_download
#SBATCH --output=logs/exp6_download_%j.out
#SBATCH --error=logs/exp6_download_%j.err
#SBATCH --partition=mc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:01:00

echo "ERROR: compute nodes have no internet — do not download via Slurm." >&2
echo "Run on the login node in tmux:" >&2
echo "  bash experiments/download/download_all_login.sh" >&2
exit 1
