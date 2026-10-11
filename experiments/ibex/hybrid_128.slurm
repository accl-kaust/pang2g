#!/bin/bash --login
# Hybrid dispatch (priority queue + chunked pops, G2G_CHUNK) at 16-128 threads, whole Rome node.
#     sbatch experiments/ibex/hybrid_128.slurm
#SBATCH --job-name=g2g_ibex_hybrid
#SBATCH --partition=batch
#SBATCH --constraint=rome
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --output=logs/ibex_hybrid_%j.out
set -uo pipefail
source "$SLURM_SUBMIT_DIR/experiments/ibex/common.sh"
export OMP_PROC_BIND=true OMP_PLACES=cores
for C in 1 4 8 32; do
  echo "[ibex] === G2G_CHUNK=$C ==="
  G2G_CHUNK=$C THREADS="1 16 32 64 128" REPS=3 bash demo/07_parallel_scaling.sh 2>&1 | tee logs/07_chunk$C.log | tail -12
  cp demo/work/parallel/parallel_scaling.tsv results/ibex/hybrid_chunk${C}.tsv
done
mkdir -p "$SLURM_SUBMIT_DIR/results/ibex" && cp results/ibex/hybrid_chunk*.tsv logs/07_chunk*.log "$SLURM_SUBMIT_DIR/results/ibex/"
echo "[ibex] $(date) DONE_HYBRID"
