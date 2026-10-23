# G2G Phase 3 — C++ Core Implementation Plan

> **Goal**: A correct, single-threaded C++ system that produces the same
> output as the Python prototype on all test cases.
> **Duration**: Weeks 7–18 (approx. 10–12 weeks).
> **Entry state**: Phase 2 complete (profiling data in `results/`).
> **Phase 4 hook**: `GraphAligner::align()` is the entry point Phase 4 will
> parallelise; keep the snarl loop isolated in `align_snarls()`.

---

## Quick orientation: what exists vs what needs writing

| File | Status | Notes |
|------|--------|-------|
| `include/g2g/types.hpp` | ✅ Ready | All core structs defined |
| `include/g2g/graph.hpp` | ✅ Ready | CSR layout specified |
| `src/graph.cpp` | ✅ Ready | `add_node`, `add_path`, `finalise` (CSR builder) — **fully implemented** |
| `include/g2g/sequence.hpp` | ✅ Ready | |
| `src/sequence.cpp` | ✅ Ready | `edit_distance` (Wagner-Fischer) — **fully implemented** |
| `include/g2g/matching.hpp` | ✅ Ready | |
| `src/matching.cpp` | ⚠️ Partial | `greedy_match` works; `brute_match` enum incomplete; `hungarian_match` = greedy stub |
| `include/g2g/decompose.hpp` | ✅ Ready | vEB layout design documented |
| `src/decompose.cpp` | 🔴 Stub | |
| `include/g2g/anchor.hpp` | ✅ Ready | Minimizer index design documented |
| `src/anchor.cpp` | 🔴 Stub | |
| `include/g2g/align.hpp` | ✅ Ready | |
| `src/align.cpp` | 🔴 Stub | |
| `include/g2g/io.hpp` | ✅ Ready | |
| `src/io.cpp` | 🔴 Stub | |
| `tests/test_sequence.cpp` | ✅ Ready | All test cases written |
| `tests/test_matching.cpp` | ✅ Ready | Oracle tests written |
| `tests/test_align.cpp` | ✅ Ready | Integration test written (TODOs marked) |

---

## Task list in priority order

### Week 1 (highest priority per Phase 2)

**Task P3-1: Complete `brute_match` enumeration**
- File: `src/matching.cpp`, function `brute_match()`
- Implement the combinations/permutations loop — mirrors Python's
  `itertools.combinations` + `itertools.permutations` exactly.
- Gate on `k ≤ 6` with a runtime assert.
- Run `test_matching` and confirm `HungarianMatchesBrute_k2` passes.
- **Why first**: `brute_match` is the correctness oracle for all subsequent
  matching work. Nothing else can be verified without it.

**Task P3-2: Vendor LAPJV for `hungarian_match`**
- Copy `https://github.com/dkurt/lapjv` (MIT, ~200 lines of C++) into
  `src/lapjv.cpp` + `include/g2g/lapjv.hpp`.
- Wire into `hungarian_match()`: build cost matrix `C[i][j] = 1 - sim[i][j]`,
  pad rectangular case, call `lapjv()`, extract assignment.
- Run `test_matching` on all k=2..6 cases; verify score matches `brute_match`
  to within 1e-4.
- **Why urgent**: Phase 2 shows k≥7 takes 170ms+ per snarl pair with brute
  force; Hungarian drops this to ~0.3ms. Without it, the test suite is
  capped at tiny graphs.

### Weeks 2–3

**Task P3-3: Implement `decompose.cpp` (built-in DAG back-end)**
- Port `SnarlDecomposer` from `g2g/decompose.py` to C++.
- Fix the topological_sort call count: call it **once** per graph and cache
  in `TopoCache`. (Phase 2: `topological_sort` called 1.8M times = 6% of runtime.)
- Use NetworkX-free BFS reachability (raw adjacency list).
- Run `test_decompose` on the demo graph; verify snarl count = 2.

**Task P3-4: Implement vEB memory layout in `SnarlDecomposer::apply_veb_order()`**
- After DFS decomposition, re-sort `SnarlTree::snarls` into van Emde Boas
  recursive order.
- Algorithm: recursive layout of tree into array, split at sqrt(N) height:
  ```
  veb_layout(tree, node):
    if leaf: emit node
    else:    veb_layout(top_half)
             emit node
             for each child: veb_layout(child)
  ```
- Verify: for a balanced tree of 1000 snarls, sequential bottom-up traversal
  in vEB order should produce ≥2× fewer L3 cache misses than BFS order.
  (Measure with `perf stat -e cache-misses` on the inner DP loop microbenchmark
  from Phase 2 — this is the one remaining Phase 2 item.)

**Task P3-5: Implement `find_anchors` (naïve path first, then minimizer index)**
- Phase 3a (naïve, 1 day): direct port of `_find_anchors()` from `align.py`.
  Correct but O(|VA|·|VB|). Unblocks test_align integration tests.
- Phase 3b (minimizer index, 1 week): implement `MinimizerIndex`:
  - `build()`: for each backbone node, compute w-mer minimizers, insert into
    hash map.
  - `query()`: hash the query sequence's minimizers, intersect with index.
  - Switch `find_anchors()` to use the index when |backbone_a|·|backbone_b| > 10^8.
  - Measure speedup and anchor recall vs naïve on synthetic graphs.

### Weeks 4–6

**Task P3-6: Implement `GraphAligner::align_snarls()` and `score_snarl_pair()`**
- Direct port of `_align_snarls()` + `_score_snarl_pair()` from `align.py`.
- Use `allele_sequences()` from `sequence.cpp` to get inner node sequences.
- Call `allele_match()` (dispatcher) for the W_struct term.
- Run `test_align::AlignsTwoSnarlPairs` and verify result.

**Task P3-7: Implement `build_deltas()`**
- Fill `AlignmentResult::a_only_alleles` and `b_only_alleles`.
- Complete the TODO markers in `test_align::CTTAGGAATCGisAOnly` and
  `GATTACisBOnly`.

**Task P3-8: GFA reader**
- Parse S (segment), L (link), P (path) lines from GFA v1.1.
- Store sequences in `VariationGraph::seq_data` (flat array, no per-node string).
- Call `vg.finalise()` before returning.
- Test: round-trip a small hand-written GFA file; verify node count, edge count,
  path sequences match the original.

**Task P3-9: Output writer**
- Implement `write_alignment_gaf()` and `write_delta_vcf()`.
- Test: align demo graphs, write output, parse it back, verify round-trip.

### Weeks 7–8 (ongoing)

**Task P3-10: Test suite expansion and fuzz testing**
- For every test case in `g2g/demo.py`, compare C++ vs Python output.
- Fuzz test: generate random variation graphs, run aligner, check:
  - `total_score ≥ 0`
  - No duplicate node matches
  - `delta ∪ shared = full node set` of each graph
- Script: `scripts/oracle_check.sh` runs both implementations and diffs
  the output files.

---

## Performance targets after Phase 3

At the end of Phase 3 (before Phase 4 parallelism), the single-threaded
C++ system should achieve on the demo graph:
- Runtime: < 1ms (Python: ~75ms at n=20 — target 100× speedup).
- On a 10K-snarl synthetic graph: < 1s (Python: ~35s — target 35×).

These speedups come from:
1. Removing brute-force matching (P3-2): 1,491× at k=7.
2. SIMD edit_distance (P3-TODO): 8–16×.
3. Cached topological sort (P3-3): removes 1.8M redundant calls.
4. Minimizer index (P3-5): replaces O(|VA|·|VB|) pairwise scan.

---

## Build instructions

```bash
# Prerequisites (conda)
conda create -n g2g_cpp -c conda-forge cmake cxx-compiler
conda activate g2g_cpp

# Build
cd G2G_cpp
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j$(nproc)

# Run tests
cd build && ctest --output-on-failure

# Run demo (once GFA reader is implemented)
./build/g2g demo_a.gfa demo_b.gfa
```

---

## Phase 4 preparation (keep in mind during Phase 3)

Phase 4 (work-stealing parallel runtime) will replace the snarl loop in
`GraphAligner::align_snarls()` with a task DAG. To make this easy:

1. **Do not put side effects in `score_snarl_pair()`** — it must be a pure
   function of `(ga, gb, sa, sb)` with no shared mutable state.
2. **Annotate subtree work in `SnarlTree`** — `Snarl::subtree_work` is the
   field Phase 4 uses for the steal-largest priority function.
3. **NUMA-aware allocation**: Phase 4 will call `allele_sequences()` per-task
   and needs to allocate the result on the NUMA node of the executing thread.
   Keep `allele_sequences()` allocation-free by writing into a caller-supplied
   `std::vector<std::string>& out` buffer (already done in the header).
