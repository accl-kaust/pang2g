# g2g: the Python prototype

This is the prototype implementation of graph-to-graph alignment for
pangenome variation graphs, written before the C++ engine and kept as a
readable reference and correctness oracle. Where current tools map a
sequence (a read or an assembly) to a graph (a pangenome reference), `g2g`
maps a graph to another graph, comparing two population pangenomes directly
without routing everything through a single linear reference.

## Background

Genome reference models have gone through three stages: a single linear
string (GRCh38), where sequences are compared to a sequence; the pangenome
variation graph (vg, minigraph), where sequences are mapped to a graph; and
what this project does, aligning two pangenome graphs to each other.

Once you have population-specific pangenomes, say
one built from African cohorts and one from East Asian cohorts, you want to
ask which structural variants are shared, which are population-specific,
and where the graphs diverge. A single linear reference cannot answer those
questions without reintroducing the very reference bias that pangenomes were
built to eliminate.

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

The score is `W(A) = W_seq + W_topo + W_struct`, where `W_seq` is the
sequence similarity of matched boundary (anchor) nodes, `W_topo` the
structural class consistency (both simple bubbles, nested snarls, and so
on), and `W_struct` the bipartite matching score across the allele sets
within each snarl pair. The complexity is `O(|T| · |T'| · k³)`, where `T`
and `T'` are the snarl trees and `k` the maximum number of alleles per
snarl, typically 2 to 5 in real pangenomes.

## Installation

Python 3.12 or newer and `networkx`, nothing else. With miniconda:

```bash
conda create -n g2g python=3.12 networkx
conda activate g2g
```

## Quick start

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

The split between memory and compute files is deliberate. The memory files
(`types.py`, `graph.py`) can be read to understand the data model without
knowing anything about the algorithms, and the compute files
(`sequence.py`, `decompose.py`, `align.py`) hold all the algorithmic logic
and can be swapped independently, for instance replacing `sequence.py` with
a C extension for large graphs without touching anything else.

## Reading the output

```
ANCHOR NODES  (shared backbone)
  b1a ['ACGT']  ↔  b1a ['ACGT']   sim=1.00
```
These are nodes whose sequences are identical, or nearly so, in both graphs.
They act as alignment anchors; everything else is interpreted relative to
them.

```
SNARL ALIGNMENTS
  Snarl (s0,m1)  ↔  Snarl (s0,m1)   score=5.500
    A-allele[0] 'ACGT'  ↔  B-allele[0] 'ACGT'  sim=1.00   ← shared
    A-allele[1] 'ATGT'  ↔  B-allele[1] 'ATGT'  sim=1.00   ← shared
```
One record per aligned snarl pair. Matched alleles are homologous across
populations, and the score decomposes into `W_seq + W_topo + W_struct`.

```
DELTA GRAPH  ΔA→B  (A-only alleles)
  snarl (m1,s_end)  seq='CCCGCTTAGGAATCGTTAAG'

DELTA GRAPH  ΔB→A  (B-only alleles)
  snarl (m1,s_end)  seq='CCCGGATTACTTAAG'
```
The delta graphs are the main output. They list every allele present in one
population graph but absent from the other, the population-specific
structural variants, the graph-level analogue of a VCF.

```
TOTAL ALIGNMENT SCORE:  13.0000
```
A single number summarising overall graph similarity, higher meaning more
structurally similar. It can be used to build a distance matrix across many
population pangenomes.

## Limitations of the prototype

It handles DAGs only; graphs with cycles (tandem repeats, for example) are
not supported, and the snarl decomposer raises `ValueError` on a cyclic
graph. It only finds simple bubbles, where all paths from the source
converge at a single sink; nested structures are detected but not
recursively scored in the DP. The bipartite allele matching is brute force,
exact but exponential in `k`, so it is safe up to about `k = 8`; beyond that
one would use the Hungarian algorithm
(`scipy.optimize.linear_sum_assignment`). Snarl pairing is greedy, one pair
at a time in topological order; a global DP over the whole snarl tree would
give optimal pairings at a higher cost. And there is no GFA I/O: graphs are
built programmatically. All of these are addressed in the C++ engine in
`src/`, which reads GFA, uses LAPJV for the matching, offers interval
pairing, and runs in parallel.

## Theoretical basis

The alignment score `W(A) = W_seq + W_topo + W_struct` generalises pairwise
sequence alignment to graphs. When both graphs are paths it reduces to
standard sequence alignment (Needleman-Wunsch); when one graph is a path it
reduces to sequence-to-graph alignment (as in GraphAligner); when both are
arbitrary graphs it is this project.

The snarl tree decomposition reduces the NP-hard general graph alignment
problem to a polynomial-time tree alignment problem, exploiting the low
treewidth of real variation graphs (typically 3 to 5 or less for
bubble-decomposed human pangenomes). The docstrings in `align.py` and
`decompose.py` carry the full derivation.
