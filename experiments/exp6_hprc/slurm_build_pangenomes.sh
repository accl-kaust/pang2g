#!/bin/bash --login
# SLURM - Build HPRC pangenomes (Stage 2 of Exp6)
# Balanced panel is small (~5 samples/pop = ~6 genomes incl. GRCh38), so each
# minigraph-cactus build is a few hours; both fit inside the cluster's 14h
# walltime cap. If a build ever bumps the limit, split AFR/EAS into two jobs.
#SBATCH --job-name=g2g_exp6_build
#SBATCH --output=logs/exp6_build_%j.out
#SBATCH --error=logs/exp6_build_%j.err
#SBATCH --partition=batch      # adjust to your cluster
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
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

echo "[slurm] Build job $SLURM_JOB_ID started at $(date)"

# Re-stratify first so the sample lists contain ONLY samples whose assemblies
# are actually on disk (an assembly download can lag/fail), and stay balanced
# AFR/EAS. This rewrites *_samples.txt before the build reads them.
echo "[slurm] Stratifying (existence-filtered, balanced)..."
python3 experiments/exp6_hprc/01_stratify_hprc_samples.py \
    --panel        "$G2G_DATA_ROOT/hprc/metadata/1kgp_pops.panel" \
    --index        "$G2G_DATA_ROOT/hprc/metadata/hprc_v2_index.tsv" \
    --out-dir      "$G2G_DATA_ROOT/hprc/metadata/" \
    --assembly-dir "$G2G_DATA_ROOT/hprc/assemblies/"

bash experiments/exp6_hprc/02_build_pangenomes.sh
echo "[slurm] Build completed at $(date)"
