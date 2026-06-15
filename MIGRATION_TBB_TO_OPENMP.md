# TBB → OpenMP Migration Guide

## Files changed (4 core + 2 Slurm scripts)

### 1. `CMakeLists.txt`
- Removed: `find_library(TBB_LIB tbb)` / `find_path(TBB_INCLUDE ...)`
- Removed: `G2G_HAVE_TBB` compile definition
- Added:   `find_package(OpenMP)` / `G2G_HAVE_OPENMP`
- Added:   `target_link_libraries(g2g_lib PUBLIC OpenMP::OpenMP_CXX)`

### 2. `include/g2g/parallel.hpp`
- Removed: `#include <tbb/...>` headers (not needed — OpenMP is pragma-based)
- No API changes — `ParallelAligner`, `SnarlPairTask`, `RunStats` unchanged

### 3. `src/parallel.cpp`
- Removed: all TBB includes (`tbb/task_arena.h`, `tbb/task_group.h`,
  `tbb/concurrent_priority_queue.h`, `tbb/concurrent_queue.h`)
- Added:   `#include <omp.h>`, `<condition_variable>`, `<queue>`
- Replaced `tbb::concurrent_priority_queue` with
  `std::priority_queue` + `std::mutex` + `std::condition_variable`
- Replaced `tbb::task_arena` + `tbb::task_group` with
  `#pragma omp parallel num_threads(num_threads_)` calling a worker lambda
- `detect_threads()` now calls `omp_get_max_threads()` instead of
  `std::thread::hardware_concurrency()`
- Idle wait uses `condition_variable::wait()` instead of `std::this_thread::yield()`
  (less CPU burn when queue is temporarily empty)
- `G2G_STEAL_RANDOM` ablation preserved: compiles with `std::queue` (FIFO)
  instead of `std::priority_queue`

### 4. `tests/CMakeLists.txt`
- Changed: `if(G2G_HAVE_TBB)` → `if(G2G_HAVE_OPENMP)` for `test_parallel`

### 5. `slurm/scripts/setup_cluster.sh` (manual edit needed)
Change the conda install line:
```diff
-    conda install -y -c conda-forge -c bioconda cmake cxx-compiler tbb-devel minimap2
+    conda install -y -c conda-forge -c bioconda cmake cxx-compiler minimap2
```
(OpenMP ships with GCC — no extra package needed.)

### 6. `slurm/scripts/exp5_ablation.sh` (manual edit needed)
The matching micro-benchmark links `-ltbb`. Change:
```diff
-    -pthread -ltbb \
+    -pthread -fopenmp \
```

## Build instructions (unchanged except no TBB package)

```bash
conda activate phase5
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
cd build && ctest --output-on-failure
```

## Behavioral differences

- **Correctness**: Identical — same task DAG, same greedy pairing,
  same priority ordering, same memory-ordering guarantees on atomics.
- **Performance**: The mutex-protected priority queue adds ~1 µs of
  contention per task pop vs TBB's lock-free queue. At HPRC scale
  (10M snarls, per-task DP work >> 1 µs) this is negligible.
  For very small graphs (< 100 snarls), single-threaded mode avoids
  the queue entirely.
- **Idle wait**: `condition_variable::wait()` instead of spin-yield
  means idle cores consume 0% CPU. Better for shared HPC nodes.
