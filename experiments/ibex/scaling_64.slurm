#!/bin/bash --login
# Parallel scaling to 128 threads on a whole IBEX Rome node (EPYC 7702), plus the ablation.
#     sbatch experiments/ibex/scaling_128.slurm
#SBATCH --job-name=g2g_ibex_scale64
#SBATCH --partition=batch
#SBATCH --constraint=rome
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --output=logs/ibex_scale64_%j.out
set -uo pipefail
source "$SLURM_SUBMIT_DIR/experiments/ibex/common.sh"
export OMP_PROC_BIND=true OMP_PLACES=cores
THREADS="1 2 4 8 16 32 64" REPS=3 bash demo/07_parallel_scaling.sh 2>&1 | tee logs/07_ibex.log | tail -14
cp demo/work/parallel/parallel_scaling.tsv results/ibex/parallel_scaling_64.tsv
PTHREADS=32 SCHED_THREADS=8 bash demo/08_ablation.sh 2>&1 | tee logs/08_ibex.log | tail -12
cp demo/work/ablation/ablation_results.tsv results/ibex/ablation_results_64.tsv
mkdir -p "$SLURM_SUBMIT_DIR/results/ibex" && cp results/ibex/* logs/*.log "$SLURM_SUBMIT_DIR/results/ibex/"
echo "[ibex] $(date) DONE_SCALE"
