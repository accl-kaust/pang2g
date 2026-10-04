#!/bin/bash --login
# Camera-ready re-measurements on the ACCL cluster (EPYC 7763 nodes), part 1:
# build + tests, pairing-error sweep, parallel scaling (1-48 threads: the build
# partition's limit) and the ablation.  Submit from the repo root:
#     sbatch experiments/accl/cr_build_measure.slurm
#SBATCH --job-name=g2g_cr_measure
#SBATCH --partition=build
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --time=08:00:00
#SBATCH --output=logs/cr_measure_%j.out
#SBATCH --error=logs/cr_measure_%j.err
set -uo pipefail
export PATH=$HOME/miniconda3/envs/g2g/bin:$HOME/miniconda3/bin:/opt/slurm/current/bin:$PATH
export OMP_PROC_BIND=true OMP_PLACES=cores
R="${SLURM_SUBMIT_DIR:-$PWD}"; cd "$R"; mkdir -p logs results/accl
echo "[cr] $(date) host=$(hostname) cpus=$SLURM_CPUS_ON_NODE mem=$(free -g | awk '/Mem/{print $2}')G gcc=$(gcc --version | head -1)"

echo "[cr] === build ==="
# cmake/ninja are not installed on every node; fall back to the copy Vitis
# bundles (/tools/Xilinx is mounted on all nodes) and to Unix Makefiles.
command -v cmake >/dev/null || export PATH="/tools/Xilinx/Vitis/2022.2/tps/lnx64/cmake-3.21.4/bin:$PATH"
command -v cmake >/dev/null || { echo "[cr] no cmake on $(hostname)"; exit 1; }
if command -v ninja >/dev/null; then GEN=Ninja; else GEN="Unix Makefiles"; fi
echo "[cr] cmake=$(command -v cmake) ($(cmake --version | head -1)) generator=$GEN"
rm -rf build/release
cmake -B build/release -G "$GEN" -DCMAKE_BUILD_TYPE=Release . > logs/cmake_configure.log 2>&1 || { tail -20 logs/cmake_configure.log; exit 1; }
cmake --build build/release -j"$SLURM_CPUS_ON_NODE" > logs/cmake_build.log 2>&1 || { grep -n -A5 "error" logs/cmake_build.log | head -40; exit 1; }
ln -sfn release/g2g build/g2g
(cd build/release && ctest -j8 2>&1 | tail -4)

echo "[cr] === pairing-error sweep (exp1b) ==="
python3 experiments/exp1_correctness/pairing_error_sweep.py --g2g build/release/g2g \
    --pairing both --anchor-lens 8,2,1 --private-fracs 0,0.1 --linker-every 4 \
    --out results/accl/exp1b_pairing_error.csv > logs/exp1b.log 2>&1 && tail -2 logs/exp1b.log

echo "[cr] === parallel scaling (demo/07, 1-48 threads, 3 reps) ==="
THREADS="1 2 4 8 16 32 48" REPS=3 bash demo/07_parallel_scaling.sh 2>&1 | tee logs/07_parallel_scaling_accl.log | tail -14
cp demo/work/parallel/parallel_scaling.tsv results/accl/parallel_scaling_48.tsv

echo "[cr] === single-thread phase breakdown on the scaling pair (new decomposer) ==="
for i in 1 2 3; do build/release/g2g demo/work/parallel/demoA.gfa demo/work/parallel/demoB.gfa --threads 1 --verbose-timing --gaf /dev/null 2>&1 | grep -E "timing|runstats"; done | tee results/accl/phase_breakdown_1thread.txt

echo "[cr] === ablation (demo/08, PTHREADS=32, SCHED_THREADS=8) ==="
PTHREADS=32 SCHED_THREADS=8 bash demo/08_ablation.sh 2>&1 | tee logs/08_ablation_accl.log | tail -12
cp demo/work/ablation/ablation_results.tsv results/accl/ablation_results.tsv

echo "[cr] $(date) DONE_MEASURE"
