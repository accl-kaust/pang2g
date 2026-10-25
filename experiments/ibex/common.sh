# Sourced by the IBEX jobs: modules, per-job copy of the tree, native build, tests.
module load gcc/12.2.0 cmake/3.28.4/gnu-11.3.1 python/3.11.0 2>&1 | tail -1
export PATH=/ibex/user/$USER/venv/bin:$PATH
SRC="${SLURM_SUBMIT_DIR:?}"
JOB="/ibex/user/$USER/run_${SLURM_JOB_ID}"
mkdir -p "$JOB" && rsync -a --exclude build --exclude demo/work --exclude .git "$SRC/" "$JOB/" && cd "$JOB"
mkdir -p logs results/ibex
echo "[ibex] $(date) host=$(hostname) $(lscpu | grep 'Model name' | sed 's/  */ /g') cpus=$(nproc) mem=$(free -g | awk '/Mem/{print $2}')G"
GT=/ibex/user/$USER/googletest-1.14.0; GTARG=""; [ -d "$GT" ] && GTARG="-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=$GT"
cmake -B build/release -DCMAKE_BUILD_TYPE=Release $GTARG . > logs/cmake_configure.log 2>&1 || { tail -20 logs/cmake_configure.log; exit 1; }
cmake --build build/release -j"$(nproc)" > logs/cmake_build.log 2>&1 || { grep -n -A5 "error" logs/cmake_build.log | head -40; exit 1; }
ln -sfn release/g2g build/g2g
(cd build/release && ctest -j16 2>&1 | tail -4)
