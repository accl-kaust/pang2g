#!/bin/bash --login
# Reference-anchored single-core sweep from 4 to 8 Gbp on a 730 GB Turin node.
#     sbatch experiments/ibex/scaling_8gbp.slurm
#SBATCH --job-name=g2g_ibex_8gbp
#SBATCH --partition=batch
#SBATCH --constraint=turin
#SBATCH --cpus-per-task=8
#SBATCH --mem=700G
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --output=logs/ibex_8gbp_%j.out
set -uo pipefail
source "$SLURM_SUBMIT_DIR/experiments/ibex/common.sh"
START_MB=4096 MAX_MB=8192 MEM_CEIL_MB=650000 BUDGET_SEC=40000 bash demo/06_scaling_limit.sh 2>&1 | tee logs/06_ibex.log | tail -20
mkdir -p "$SLURM_SUBMIT_DIR/results/ibex"
cp demo/work/scaling_limit/scale_limit_results.tsv "$SLURM_SUBMIT_DIR/results/ibex/scale_limit_results_8gb.tsv"; cp logs/06_ibex.log "$SLURM_SUBMIT_DIR/results/ibex/"
echo "[ibex] $(date) DONE_8GBP"
