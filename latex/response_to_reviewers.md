# Response to the reviewers — BIBM 2026 paper "PanG2G: Scalable Base-Pair-Resolution Graph-to-Graph Alignment of Pangenome Graphs"

We thank the three reviewers for their careful reading. Below we answer every
point and say where the camera-ready changed. Section and figure numbers refer
to the camera-ready version. Beyond the requested edits, the revision
re-measured two claims of the submitted version and corrected a decomposer
bug that the new real-data check exposed; both are stated openly in the paper
(Q1) and summarized at the end of this letter.

*(The reviewer comments are paraphrased from our notes; please substitute the
original wording before sending.)*

---

## Reviewer 1

**R1.1 — The problem is never formally defined.**
Added §III-A "Problem Definition" with two definitions: variation graph /
snarl / allele (Def. 1) and graph-to-graph alignment as the maximization of
Eq. (1) over node, edge and snarl maps that preserve snarl-tree ancestry and
backbone order (Def. 2). The text now states what PanG2G solves exactly (the
within-snarl allele matching) and what it solves heuristically (the snarl map),
and that the delta graph is the set of alleles left unmatched, including every
allele of a snarl that has no partner.

**R1.2 — Related work is a tool catalogue; position the problem against
graph alignment more broadly.**
§II was rewritten around three questions: global alignment beyond strings
(Needleman-Wunsch, partial-order alignment of two POGs [Lee et al. 2002], the
NP-completeness of aligning alignments [Kececioglu & Starrett 2004] and of
sequence-to-graph alignment with graph edits [Jain et al. 2020]); general graph
comparison (graph edit distance via bipartite assignment [Riesen & Bunke 2009],
network alignment [IsoRank]); and pangenome graph comparison (vg deconstruct as
used by the HPRC, ODGI path similarity, PanGraph, SplitMEM, PGVF). Read mappers
appear only as the linearize-and-map baseline. Six references were added and
checked (DOIs, author lists; the ParSGA entry now carries its full author list).

**R1.3 — Code and data availability.**
The availability paragraph (§VI) now points to the public repository, which
contains the C++ engine, the Python prototype, the synthetic generators with
seeds, every experiment harness, the measured result tables behind every
figure, and a single `make_figures.py` that regenerates all figures from those
tables. It names every public input with its version (1000 Genomes phase 3
release 20130502, GRCh37, gnomAD-SV v2.1, ClinVar GRCh37 weekly release of
28 Sep 2026, PGGB HLA-DRB1 test haplotypes). The repository carries an MIT
license, a README with build/run instructions, and a `CITATION.cff`; the
camera-ready release is archived on Zenodo (DOI given in §VI).

## Reviewer 2

**R2.1 — No code or data was available to assess the claims.**
See R1.3. Everything reported, including the new experiments below, is
reproducible from the repository; `experiments/` and `demo/` contain the
harnesses, and `results/` the measured tables.

**R2.2 — Greedy snarl pairing: a ~20% error on toy graphs is high, and errors
may compound with graph size; the paper should measure error as a function of
size.**
We agree, and this point led to the largest change. Three things happened:

1. *The "~20%" figure was not a measurement.* It was the 0.8 assertion
   threshold of our unit tests. Re-measuring on the planted graphs of the
   submission (8-bp boundary nodes) the engine recovers the planted deltas with
   precision 1.000 and recall 0.999–1.000 at every size from 5 to 10^5 snarls
   (Q1, "Planted graphs"). The camera-ready says so explicitly.
2. *The real failure mode was identified and measured.* The default pairing
   rule matches snarls by identical boundary (source, sink) sequences and
   resolves repeated keys positionally. It fails exactly when keys repeat and a
   snarl exists in only one graph, which shifts every later snarl with that
   key. On real vg construct graphs 10–13% of bubbles share a key with another
   bubble (35–40% of boundary nodes are ≤2 bp), so this matters. Q1, "Pairing
   error at scale" and the new Fig. 4 measure it on planted graphs with 1- and
   2-bp boundary nodes and 10% private snarls, from 5 to 10^5 snarls: the
   delta F1 of the boundary rule falls from 1.00 to 0.55 (precision 0.38) at
   10^5 snarls, i.e. the compounding the reviewer suspected is real.
3. *The proposed fix was implemented and evaluated.* `--pairing interval`
   chains the shared anchors in topological order and pairs the snarls between
   consecutive anchors with an order-preserving dynamic program (§III-B). It
   holds F1 ≥ 0.999 in every configuration up to 10^5 snarls at no measurable
   cost, and on the real HLA-DRB1 graphs pairs 972 snarls instead of 694,
   leaving 192+52 instead of 470+330 unpaired. Unpaired snarls now contribute
   all their alleles to the delta (as Def. 2 requires; the submitted code
   silently dropped them).

The default remains the boundary rule so that the scaling and ablation numbers
are reproduced as reported; the Limitations state that the interval rule is
optimal only within an interval and that full tree matching remains the
principled end point.

**R2.3 — Why a single centralized queue rather than work stealing? Explain
the 128-thread efficiency loss.**
§III-D ("Centralized priority-dispatch runtime") now gives the reason: a single
LPT queue preserves largest-first ordering, which randomized stealing loses, at
the price of a lock that may limit scaling at high thread counts. The
camera-ready restricts the scaling evaluation and its claims to 64 cores and
re-measures it there (51.9x DP-phase speedup at 64 cores, 81% efficiency,
Fig. 5a); it makes no claim about higher thread counts. We implemented the
hybrid we proposed (a worker takes a chunk of the largest ready tasks per lock
acquisition) and measured it at 64 cores: chunks of 4, 8 and 32 give 53.4x,
48.8x and 15.2x against 52.0x for the plain queue, so it brings no gain there
(Q2). The ablations show that ordering (FIFO 1.20x) and dynamic dispatch
(static 1.92x) both matter.

## Reviewer 3

**R3.1 (major) — Variant representation differs across builders
(Minigraph-Cactus vs. PGGB); coordinate keying ignores normalization.**
§IV-C now explains why same-builder comparison is a precondition: beyond
indel left-alignment, builders cut regions into different snarls, nest
variation differently, merge/split adjacent SNPs and place repeat variants
differently; neither mode reconciles this, and all comparisons in the paper
are same-builder. The new HLA experiment (Q1, "Real vg construct graphs")
quantifies the effect even within one builder: recall 0.99 but precision 0.78,
every false positive being a differently placed indel/SNP cluster.

**R3.2 (major) — The anchor definition (degree ≤ 1 nodes) is unclear; how
does node chopping affect it?**
§III-B "Anchor detection" now defines anchors precisely and reports what the
rule does on real vg construct graphs: half the nodes satisfy the degree rule,
but over 80% of those are 1-bp SNP allele nodes (the rule admits allele nodes
of simple bubbles, so shared alleles act as anchors too); nodes shorter than
the 15-bp minimizer length are never anchored, so anchors come from the ≤32-bp
reference chunks between variants. We also state honestly that anchors
contribute only the α term in the default pairing and pin the coordinate
systems only in the new anchor-interval pairing; the submitted text implied
more.

**R3.3 (major) — Depth-1 scope: are nested Minigraph-Cactus graphs handled?**
No. §III-B states that nested snarls are detected but only top-level snarls
are kept (depth-1), that no Minigraph-Cactus results are reported, and that
the nested recurrence, vEB layout and DAG dependencies are exercised only at
depth 1; the scripts for HPRC/Minigraph-Cactus graphs are in the repository
as scripts only.

**R3.4 (major) — Built-in snarl finder vs. `vg snarls`: what is the
difference, and does the vg path work?**
Both backends are described in §III-B and both were run on the real HLA graphs
(Q1). This check exposed a bug in the submitted finder: its "first node reached
by all successors" rule is not a separator test and closed a bubble one node
early wherever vg construct links the alternate node of one SNP directly to the
reference node of the next, dropping the following variant (recall 0.74 on the
HLA graphs). The camera-ready replaces it with a superbubble scan after Onodera
et al. (2013); the built-in and `vg snarls` backends now produce byte-identical
results (recall 0.99), the built-in finder is ~20× faster than before, and the
synthetic results are unchanged. The Q4 chromosome-22 run was repeated with the
corrected finder; its numbers changed slightly and are updated in the paper.

**R3.5 (major) — Compare against a semi-linear baseline (e.g. Giraffe).**
Q2 now discusses it: mapping each haplotype of one population to the other
population's graph removes reference bias on the mapping side but still expands
one graph into h sequences, so its cost should also grow linearly in h and it
still yields per-haplotype alignments that need calling and merging. We have
not benchmarked it and say so.

**R3.6 (major) — 20% pairing error is high.**
See R2.2: the 20% was a test threshold; the measured error is ≤0.1% on the
submitted graphs, the real failure mode is measured up to 10^5 snarls, and the
anchor-interval pairing removes it.

**R3 minor points.**
- W_seq, W_topo, the unmatched sets and the shared alignment are defined in the
  caption and body of Algorithm 1; the "shared graph" is now consistently
  called the shared alignment and described as a GAF-style file, not a
  materialized graph (abstract, Fig. 1, Algorithm 1).
- k is defined at first use (Def. 1 and the introduction).
- "Hungarian" was replaced by LAPJV (Jonker-Volgenant) throughout; the
  introduction explains what LAPJV computes on first mention.
- "Cache-oblivious" is used consistently; the vEB figure was dropped for space,
  the design is described in §III-D.
- The Q4 takeaway was reworded: the delta co-localizes with structural and
  clinical variation in general, not preferentially with population-
  differentiated SVs, because the 1000 Genomes panel is SV-depleted.
- Typos fixed ("aligments", "Limiations", "this are", stray parentheses,
  remaining British spellings).

---

## Changes made on our own initiative

- **Corrected claim (Q1).** The submitted "≥80% precision/recall, ~20% pairing
  error" was the unit-test threshold; the measured values are given instead.
- **Decomposer bug fixed** (see R3.4); real-data runs repeated.
- **Unpaired snarls** now contribute their alleles to the delta (Def. 2).
- **k=2 fraction** in Q2 corrected from ~67% (prototype workload) to 50% (the
  scaling workload).
- **New figure** (Fig. 4, pairing error vs. size); the vEB schematic, the
  F_ST-sweep figure and the scaling table were removed and their content kept
  in the text to meet the 8-page limit; scaling and throughput share one
  figure, as do the Gini and matching-cost plots.
- **Reproducibility.** Every figure is generated from committed measurement
  tables by one script; the result tables of the repeated runs are in the
  repository.
- **Re-measurements on the original node class** (`results/accl/`): with the
  corrected decomposer the throughput comparison against linearize-and-map
  improves from 12.7× to 51× at h=200 (PanG2G 7.6–12.3 s over h=2–200), the
  serial set-up of the parallel aligner drops from 7.8 s to 3.5 s (end-to-end
  speedup at 32 threads 2.59× → 3.59×), and the ablation, pairing sweep and Q4
  numbers reproduce. The 64/128-thread and 8-Gbp points could not be repeated
  within the cluster's current 48-CPU / 192 GB job limit and are reported from
  the submitted run, as stated in the text.
