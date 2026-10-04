#!/bin/bash --login
# Camera-ready re-measurement, part 2: the reference-anchored single-core
# scaling sweep (demo/06) with the superbubble decomposer, doubling from 8 Mbp
# until the 192 GB partition limit (the 4-Gbp point needs ~121 GiB).
#     sbatch experiments/accl/cr_scaling_limit.slurm
#SBATCH --job-name=g2g_cr_scaling
#SBATCH --partition=build
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --time=12:00:00
#SBATCH --output=logs/cr_scaling_%j.out
#SBATCH --error=logs/cr_scaling_%j.err
set -uo pipefail
export PATH=$HOME/miniconda3/envs/g2g/bin:$HOME/miniconda3/bin:/opt/slurm/current/bin:$PATH
R="${SLURM_SUBMIT_DIR:-$PWD}"; cd "$R"; mkdir -p logs results/accl
echo "[cr] $(date) host=$(hostname) mem=$(free -g | awk '/Mem/{print $2}')G"
# wait for the build job to have produced the binary (part 1), up to 2 h
for i in $(seq 1 240); do [ -x build/release/g2g ] && break; sleep 30; done
[ -x build/release/g2g ] || { echo "no binary"; exit 1; }
ln -sfn release/g2g build/g2g
START_MB=8 MAX_MB=4096 MEM_CEIL_MB=180000 BUDGET_SEC=40000 bash demo/06_scaling_limit.sh 2>&1 | tee logs/06_scaling_limit_accl.log | tail -20
cp demo/work/scaling_limit/scale_limit_results.tsv results/accl/scale_limit_results_4gb.tsv
echo "[cr] $(date) DONE_SCALING"
