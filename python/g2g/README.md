# g2g - Graph-to-Graph Pangenome Alignment

A prototype implementation of graph-to-graph alignment for pangenome variation graphs.
Where current tools map a *sequence* (read or assembly) to a *graph* (pangenome reference),
`g2g` maps a *graph to another graph* - comparing two population pangenomes directly,
without routing everything through a single linear reference.

---

## Background

Genome reference models have evolved in three steps:

| Era | Paradigm | Reference structure |
|-----|----------|---------------------|
| 1 | Sequence → Sequence | Single linear string (GRCh38) |
| 2 | Sequence → Graph | Pangenome variation graph (vg, minigraph) |
| 3 | **Graph → Graph** | **Two pangenome graphs aligned to each other** ← this project |

The motivation: once you have population-specific pangenomes (e.g. one built from African
cohorts, one from East Asian cohorts), you need a way to ask - which structural variants
are shared? Which are population-specific? Where do the graphs diverge? A single linear
reference cannot answer these questions without reintroducing the very reference bias that
pangenomes were built to eliminate.

---

## How it works

The alignment proceeds in five steps:

```
Graph A ──┐                          ┌── ΔA→B  (A-only alleles)
           ├─ anchor detection ──────┤
Graph B ──┘     │                    └── ΔB→A  (B-only alleles)
                │
         snarl decomposition
         (bubble finding on both graphs)
                │
         snarl pairing
         (match bubbles between graphs)
                │
         per-snarl DP
         (bipartite allele matching)
                │
         AlignmentResult
         (score + shared graph + delta graphs)
```

**Scoring:** `W(A) = W_seq + W_topo + W_struct`

- `W_seq` - sequence similarity of matched boundary (anchor) nodes
- `W_topo` - structural class consistency (both simple bubbles, nested snarls, etc.)
- `W_struct` - bipartite matching score across allele sets within each snarl pair

**Complexity:** `O(|T| · |T'| · k³)` where `T`, `T'` are snarl trees and `k` is the
maximum number of alleles per snarl (typically 2–5 in real pangenomes).

---

## Installation

Requires Python ≥ 3.12 and `networkx`. With miniconda:

```bash
conda create -n g2g python=3.12 networkx
conda activate g2g
```

No other dependencies.

---

## Quickstart

```bash
# From the directory containing the g2g/ folder:
python -m g2g.demo
```

Or from Python:

```python
from g2g import VariationGraph, GraphAligner, print_result

ga = VariationGraph("Population-A")
ga.add_node("s0", "TTAGC")
ga.add_node("b1a", "ACGT")
ga.add_node("b1b", "ATGT")
ga.add_node("s_end", "TTAAG")
ga.add_path("hap1", ["s0", "b1a", "s_end"])
ga.add_path("hap2", ["s0", "b1b", "s_end"])

# build gb similarly ...

aligner = GraphAligner(alpha=1.0, beta=0.5, gamma=2.0, gap=0.5)
result  = aligner.align(ga, gb)
print_result(result)
```

---

## Package structure

```
g2g/
├── __init__.py      re-exports all public symbols
│
├── types.py         MEMORY - pure data structures
│   │                  No computation; only holds state.
│   ├── Snarl            (source, sink, allele_paths, children)
│   ├── SnarlAlignment   (matched snarl pair + per-component scores)
│   └── AlignmentResult  (anchors, snarl alignments, delta graphs)
│
├── graph.py         MEMORY - variation graph storage
│   │                  Wraps networkx DiGraph; builders + accessors only.
│   └── VariationGraph   add_node / add_edge / add_path / seq / describe
│
├── sequence.py      COMPUTE - string-level operations
│   │                  Stateless functions; no imports from rest of package.
│   ├── edit_distance(a, b)          Wagner-Fischer DP
│   ├── seq_similarity(a, b)         normalised similarity in [0, 1]
│   └── path_sequence(graph, path)   concatenate node seqs along a path
│
├── decompose.py     COMPUTE - snarl / bubble detection
│   │                  Reads graph topology, writes Snarl objects.
│   └── SnarlDecomposer(vg)
│           .decompose() -> list[Snarl]
│
├── align.py         COMPUTE - alignment DP
│   │                  Imports from all compute modules; writes nothing to disk.
│   ├── bipartite_allele_match(seqs_a, seqs_b, ...)
│   └── GraphAligner(alpha, beta, gamma, gap, min_anchor_sim)
│           .align(ga, gb) -> AlignmentResult
│
├── report.py        I/O - result formatting
│   └── print_result(result)
│
└── demo.py          ENTRY POINT
    ├── build_demo_graphs() -> (VariationGraph, VariationGraph)
    └── main()
```

The memory / compute split is intentional:
- **Memory files** (`types.py`, `graph.py`) can be read to understand the data model with
  no knowledge of algorithms.
- **Compute files** (`sequence.py`, `decompose.py`, `align.py`) contain all algorithmic
  logic and can be swapped independently - e.g. replace `sequence.py` with a C extension
  for large graphs without touching anything else.

---

## Output explained

```
ANCHOR NODES  (shared backbone)
  b1a ['ACGT']  ↔  b1a ['ACGT']   sim=1.00
```
Nodes whose sequences are identical (or near-identical) in both graphs.
They act as alignment anchors - everything else is interpreted relative to them.

```
SNARL ALIGNMENTS
  Snarl (s0,m1)  ↔  Snarl (s0,m1)   score=5.500
    A-allele[0] 'ACGT'  ↔  B-allele[0] 'ACGT'  sim=1.00   ← shared
    A-allele[1] 'ATGT'  ↔  B-allele[1] 'ATGT'  sim=1.00   ← shared
```
Per-snarl alignment records. Matched alleles are homologous across populations.
The score decomposes into `W_seq + W_topo + W_struct`.

```
DELTA GRAPH  ΔA→B  (A-only alleles)
  snarl (m1,s_end)  seq='CCCGCTTAGGAATCGTTAAG'

DELTA GRAPH  ΔB→A  (B-only alleles)
  snarl (m1,s_end)  seq='CCCGGATTACTTAAG'
```
The delta graphs are the main output. They list every allele present in one
population graph but absent from the other - the population-specific structural
variants. This is the graph-level analogue of a VCF file.

```
TOTAL ALIGNMENT SCORE:  13.0000
```
A single number summarising overall graph similarity. Higher = more structurally
similar. Can be used to build a distance matrix across many population pangenomes.

---

## Key limitations (prototype)

- **DAG only.** Graphs with cycles (e.g. tandem repeats) are not supported.
  The snarl decomposer will raise `ValueError` on a cyclic graph.
- **Simple bubbles only.** The decomposer finds bubbles where all paths from
  source converge at a single sink. Complex nested structures are detected but
  not recursively scored in the DP.
- **Brute-force allele matching.** The bipartite matching is exact but
  exponential in `k`. Safe for `k ≤ ~8`; beyond that, replace with the
  Hungarian algorithm (`scipy.optimize.linear_sum_assignment`).
- **Greedy snarl pairing.** Snarls are paired one-by-one in topological order.
  A global DP over the full snarl tree would give optimal pairings at higher cost.
- **No GFA I/O.** Graphs are built programmatically. A `VariationGraph.from_gfa()`
  reader is the natural next step to connect to real pangenome data.

---

## Roadmap

| Priority | Feature |
|----------|---------|
| High | GFA reader (`VariationGraph.from_gfa`) |
| High | Hungarian algorithm for large allele sets |
| Medium | Full snarl-tree DP (nested snarl scoring) |
| Medium | k-mer index for anchor detection (replace backbone heuristic) |
| Low | GAF-like output format for graph alignments |
| Low | Multi-graph alignment (K > 2 pangenomes simultaneously) |

---

## Theoretical basis

The alignment score `W(A) = W_seq + W_topo + W_struct` generalises pairwise
sequence alignment to graphs. Special cases:

- Both graphs are paths → reduces to standard sequence alignment (Needleman-Wunsch)
- One graph is a path → reduces to sequence-to-graph alignment (GraphAligner)
- Both graphs are arbitrary → this project

The snarl tree decomposition reduces the NP-hard general graph alignment problem
to a polynomial-time tree alignment problem, exploiting the low treewidth of
real variation graphs (typically ≤ 3–5 for bubble-decomposed human pangenomes).

For the full theoretical derivation see the inline docstrings in `align.py` and
`decompose.py`.
