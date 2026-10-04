cmake_minimum_required(VERSION 3.18)
project(g2g VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Build type defaults ────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# ── Compiler flags ────────────────────────────────────────────
add_compile_options(-Wall -Wextra -Wpedantic)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_compile_options(-O3 -march=native -DNDEBUG)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-g -O0 -fsanitize=address,undefined)
  add_link_options(-fsanitize=address,undefined)
endif()

# ── Ablation build flags ──────────────────────────────────────
# Each flag disables one optimization to isolate its contribution.  The option
# names are paper-facing; several map to a differently-named internal macro (the
# toggle each optimization actually reads), so the mapping is set explicitly.
option(G2G_BFS_LAYOUT       "Use BFS snarl-tree layout instead of van Emde Boas" OFF)
option(G2G_STEAL_RANDOM     "Use random/FIFO steal instead of largest-task-first" OFF)
option(G2G_NAIVE_ANCHOR     "Force naïve O(|VA||VB|) anchor detection"            OFF)
option(G2G_BRUTE_MATCHING   "Force brute-force bipartite matching (k<=8)"         OFF)
option(G2G_STATIC_PARTITION "Use static work partition instead of work-stealing"  OFF)

if(G2G_BFS_LAYOUT)
  add_compile_definitions(G2G_DISABLE_VEB)             # decompose.cpp
  message(STATUS "Ablation: BFS layout (-> G2G_DISABLE_VEB)")
endif()
if(G2G_STEAL_RANDOM)
  add_compile_definitions(G2G_STEAL_RANDOM)            # parallel.cpp
  message(STATUS "Ablation: random steal (-> G2G_STEAL_RANDOM)")
endif()
if(G2G_NAIVE_ANCHOR)
  add_compile_definitions(G2G_DISABLE_MINIMIZER_INDEX) # anchor.cpp
  message(STATUS "Ablation: naive anchor (-> G2G_DISABLE_MINIMIZER_INDEX)")
endif()
if(G2G_BRUTE_MATCHING)
  add_compile_definitions(G2G_BRUTE_MATCHING)          # matching.cpp
  message(STATUS "Ablation: brute matching (-> G2G_BRUTE_MATCHING)")
endif()
if(G2G_STATIC_PARTITION)
  add_compile_definitions(G2G_STATIC_PARTITION)         # parallel.cpp
  message(STATUS "Ablation: static work partition (-> G2G_STATIC_PARTITION)")
endif()

# ── Dependencies ──────────────────────────────────────────────
find_package(Threads REQUIRED)

# OpenMP (replaces oneTBB — works on any CPU vendor)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "Found OpenMP ${OpenMP_CXX_VERSION}")
  set(G2G_HAVE_OPENMP TRUE)
  add_compile_definitions(G2G_HAVE_OPENMP)
else()
  message(WARNING "OpenMP not found — parallel runtime will not be built")
  set(G2G_HAVE_OPENMP FALSE)
endif()

# Optional: detect if libhandlegraph is available
find_library(HANDLEGRAPH_LIB handlegraph HINTS ENV CONDA_PREFIX PATH_SUFFIXES lib)
find_path(HANDLEGRAPH_INCLUDE handlegraph/handle.hpp HINTS ENV CONDA_PREFIX PATH_SUFFIXES include)

if(HANDLEGRAPH_LIB AND HANDLEGRAPH_INCLUDE)
  message(STATUS "Found libhandlegraph: ${HANDLEGRAPH_LIB}")
  add_compile_definitions(G2G_HAVE_HANDLEGRAPH)
else()
  message(STATUS "libhandlegraph not found — snarl decomposition will use built-in stub")
endif()

# ── Source files ──────────────────────────────────────────────
set(G2G_SOURCES
  src/graph.cpp
  src/sequence.cpp
  src/decompose.cpp
  src/snarls_vg.cpp
  src/anchor.cpp
  src/matching.cpp
  src/align.cpp
  src/refdiff.cpp
  src/io.cpp
  src/parallel.cpp
)

add_library(g2g_lib STATIC ${G2G_SOURCES})
target_include_directories(g2g_lib PUBLIC include)

if(HANDLEGRAPH_LIB)
  target_include_directories(g2g_lib PUBLIC ${HANDLEGRAPH_INCLUDE})
  target_link_libraries(g2g_lib PUBLIC ${HANDLEGRAPH_LIB})
endif()

target_link_libraries(g2g_lib PUBLIC Threads::Threads)

if(G2G_HAVE_OPENMP)
  target_link_libraries(g2g_lib PUBLIC OpenMP::OpenMP_CXX)
endif()

# ── Executables ───────────────────────────────────────────────
add_executable(g2g src/main.cpp)
target_link_libraries(g2g PRIVATE g2g_lib)

# ── Tests ─────────────────────────────────────────────────────
# Off-by-flag so a from-source build (e.g. the demo) doesn't need to fetch
# googletest or have network access: -DG2G_BUILD_TESTS=OFF.
option(G2G_BUILD_TESTS "Build the unit tests (fetches googletest)" ON)
enable_testing()
if(G2G_BUILD_TESTS)
  add_subdirectory(tests)
endif()

# ── Tools ─────────────────────────────────────────────────────
# bench_matching: per-snarl bipartite matching cost vs. allele count k
# (LAPJV vs. brute force); writes the CSV behind fig_matching.
add_executable(bench_matching tools/bench_matching.cpp)
target_link_libraries(bench_matching PRIVATE g2g_lib)
