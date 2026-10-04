#!/bin/bash --login
# ============================================================
# SLURM — Experiment 6 (HPRC Biological Validation)
#
# This is a multi-stage pipeline. The HPRC download and
# pangenome build are separate jobs (submitted with dependencies).
#
# Usage — submit all stages in sequence:
#
#   # Stage 1: Download assemblies (~12h, network-limited)
#   JOB1=$(sbatch --parsable experiments/exp6_hprc/slurm_download.sh)
#
#   # Stage 2: Build pangenomes (~48h per pangenome)
#   JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 \
#              experiments/exp6_hprc/slurm_build.sh)
#
#   # Stage 3: Run G2G + intersect (this script)
#   JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 \
#              experiments/exp6_hprc/slurm.sh)
#
# Or run this script alone if pangenomes are already built:
#   sbatch experiments/exp6_hprc/slurm.sh
#
# ADJUST: --partition, --account for your cluster.
# ============================================================
#SBATCH --job-name=g2g_exp6_hprc
#SBATCH --output=logs/exp6_%j.out
#SBATCH --error=logs/exp6_%j.err
#SBATCH --partition=batch      # adjust to your cluster (e.g. batch on KAUST Ibex)
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

echo "[slurm] Job $SLURM_JOB_ID  CPUs=$SLURM_CPUS_ON_NODE  RAM=$(free -h | awk '/Mem/{print $2}')  date=$(date)"

# ── Stage 1: Stratify samples ─────────────────────────────────
# Classifies strictly from the authoritative 1KGP panel fetched during
# download (1kgp_pops.panel); regenerates the same balanced lists the
# download wrote, filtered to samples whose assemblies are on disk.
echo "[slurm] Stratifying HPRC samples..."
python3 experiments/exp6_hprc/01_stratify_hprc_samples.py \
    --panel        "$G2G_DATA_ROOT/hprc/metadata/1kgp_pops.panel" \
    --index        "$G2G_DATA_ROOT/hprc/metadata/hprc_v2_index.tsv" \
    --out-dir      "$G2G_DATA_ROOT/hprc/metadata/" \
    --assembly-dir "$G2G_DATA_ROOT/hprc/assemblies/"

# ── Stage 2: Run G2G on HPRC pangenomes ──────────────────────
echo "[slurm] Running G2G on HPRC pangenomes..."
bash experiments/exp6_hprc/03_run_g2g_on_hprc.sh

# ── Stage 3: Intersect with gnomAD-SV ────────────────────────
echo "[slurm] Intersecting with gnomAD-SV..."
bash experiments/exp6_hprc/04_intersect_gnomad_sv.sh

# ── Stage 4: Intersect with ClinVar (OPTIONAL) ───────────────
# ClinVar is a bonus annotation layer; the core result is the gnomAD-SV
# enrichment above. Don't let a missing/empty ClinVar BED abort the job.
CLINVAR_BED="$G2G_DATA_ROOT/clinvar/clinvar_pathogenic_svs.bed"
if [[ -s "$CLINVAR_BED" ]]; then
    echo "[slurm] Intersecting with ClinVar..."
    bash experiments/exp6_hprc/05_intersect_clinvar.sh || \
        echo "[slurm] WARNING: ClinVar step failed (non-fatal) — skipping."
else
    echo "[slurm] SKIP ClinVar: $CLINVAR_BED missing/empty (optional layer)."
fi

echo "[slurm] Exp6 completed at $(date)"
