"""
synth_generator.py - Synthetic variation graph generator for Phase 2 profiling.

Generates pairs of variation graphs with:
  - n snarls arranged in a linear chain (backbone → snarl → backbone → ...)
  - Allele counts per snarl drawn from a Zipf/power-law distribution
    (fits HPRC empirical observation: most snarls are SNPs, few are large SVs)
  - A controllable fraction of shared vs. population-specific alleles
  - Deterministic seeding for reproducibility

Usage:
    from synth_generator import generate_graph_pair
    ga, gb = generate_graph_pair(n_snarls=1000, seed=42)
"""

from __future__ import annotations
import random
import string
import sys
import os


from g2g.graph import VariationGraph


DNA = "ACGT"


def _random_seq(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(DNA) for _ in range(length))


def _mutate_seq(seq: str, error_rate: float, rng: random.Random) -> str:
    """Introduce point mutations at given rate - simulates population divergence."""
    bases = list(seq)
    for i in range(len(bases)):
        if rng.random() < error_rate:
            bases[i] = rng.choice(DNA)
    return "".join(bases)


def _zipf_k(rng: random.Random, k_min: int = 2, k_max: int = 12,
             alpha: float = 2.2) -> int:
    """
    Sample allele count k from a Zipf distribution truncated to [k_min, k_max].
    alpha=2.2 fits the HPRC human pangenome snarl-size distribution:
      ~70% snarls have k=2 (SNPs/small indels), tail reaches k~50 in SV regions.
    We cap at k_max=12 here so the brute-force prototype stays tractable.
    """
    weights = [1.0 / (i ** alpha) for i in range(k_min, k_max + 1)]
    total = sum(weights)
    r = rng.random() * total
    cumulative = 0.0
    for i, w in enumerate(weights):
        cumulative += w
        if r <= cumulative:
            return k_min + i
    return k_max


def generate_graph_pair(
    n_snarls:      int   = 100,
    shared_frac:   float = 0.7,   # fraction of alleles shared between populations
    seq_len_mean:  int   = 8,     # mean allele inner-node sequence length
    divergence:    float = 0.05,  # per-base mutation rate for non-shared alleles
    seed:          int   = 42,
) -> tuple[VariationGraph, VariationGraph]:
    """
    Build a synthetic pair of variation graphs.

    Graph topology:
        anchor_0 → [snarl_0 bubble] → anchor_1 → [snarl_1 bubble] → ... → anchor_n

    Each snarl bubble has k alleles (k drawn from a Zipf distribution).
    `shared_frac` of alleles in each snarl are identical between A and B;
    the rest are population-specific (diverged sequences).

    Returns
    -------
    (ga, gb)  two VariationGraphs with named haplotype paths
    """
    rng = random.Random(seed)

    ga = VariationGraph(f"SynthA-n{n_snarls}")
    gb = VariationGraph(f"SynthB-n{n_snarls}")

    # -- Build backbone anchors --------------------------------
    anchors: list[str] = []
    for i in range(n_snarls + 1):
        node_id = f"anc_{i}"
        seq = _random_seq(rng.randint(4, 10), rng)
        ga.add_node(node_id, seq)
        gb.add_node(node_id, seq)   # identical anchors = shared backbone
        anchors.append(node_id)

    # -- Build snarls -----------------------------------------
    # Track paths: list of (hap_a_nodes, hap_b_nodes) per snarl
    # We'll stitch them into named haplotypes at the end.
    snarl_alleles_a: list[list[str]] = []  # node IDs per allele in graph A
    snarl_alleles_b: list[list[str]] = []  # node IDs per allele in graph B

    for s_idx in range(n_snarls):
        k = _zipf_k(rng)
        n_shared = max(1, round(k * shared_frac))
        n_only_a = rng.randint(0, max(0, k - n_shared))
        n_only_b = k - n_shared - n_only_a
        if n_only_b < 0:
            n_only_b = 0

        alleles_a: list[str] = []
        alleles_b: list[str] = []
        seq_len = max(1, int(rng.gauss(seq_len_mean, seq_len_mean * 0.5)))

        # Shared alleles
        for a_idx in range(n_shared):
            seq = _random_seq(seq_len, rng)
            node_a = f"s{s_idx}_sh{a_idx}_a"
            node_b = f"s{s_idx}_sh{a_idx}_b"
            ga.add_node(node_a, seq)
            gb.add_node(node_b, seq)
            alleles_a.append(node_a)
            alleles_b.append(node_b)

        # A-only alleles
        for a_idx in range(n_only_a):
            seq = _random_seq(seq_len, rng)
            node_a = f"s{s_idx}_ao{a_idx}"
            ga.add_node(node_a, seq)
            alleles_a.append(node_a)

        # B-only alleles (diverged from a random shared allele)
        for b_idx in range(n_only_b):
            base_seq = _random_seq(seq_len, rng)
            seq = _mutate_seq(base_seq, divergence, rng)
            node_b = f"s{s_idx}_bo{b_idx}"
            gb.add_node(node_b, seq)
            alleles_b.append(node_b)

        # Ensure at least 2 alleles per side (required for a bubble)
        while len(alleles_a) < 2:
            seq = _random_seq(seq_len, rng)
            node_a = f"s{s_idx}_fill_a{len(alleles_a)}"
            ga.add_node(node_a, seq)
            alleles_a.append(node_a)
        while len(alleles_b) < 2:
            seq = _random_seq(seq_len, rng)
            node_b = f"s{s_idx}_fill_b{len(alleles_b)}"
            gb.add_node(node_b, seq)
            alleles_b.append(node_b)

        snarl_alleles_a.append(alleles_a)
        snarl_alleles_b.append(alleles_b)

    # -- Add edges: connect anchor → allele → next anchor -----
    for s_idx in range(n_snarls):
        src = anchors[s_idx]
        snk = anchors[s_idx + 1]
        for node_a in snarl_alleles_a[s_idx]:
            ga.add_edge(src, node_a)
            ga.add_edge(node_a, snk)
        for node_b in snarl_alleles_b[s_idx]:
            gb.add_edge(src, node_b)
            gb.add_edge(node_b, snk)

    # -- Register representative haplotype paths ---------------
    # One path per graph threading allele[0] of each snarl
    path_a = [anchors[0]]
    path_b = [anchors[0]]
    for s_idx in range(n_snarls):
        path_a.append(snarl_alleles_a[s_idx][0])
        path_a.append(anchors[s_idx + 1])
        path_b.append(snarl_alleles_b[s_idx][0])
        path_b.append(anchors[s_idx + 1])

    ga.paths["hap_A0"] = path_a
    gb.paths["hap_B0"] = path_b

    return ga, gb


if __name__ == "__main__":
    for n in [10, 100, 1000]:
        ga, gb = generate_graph_pair(n_snarls=n)
        print(f"n={n:5d}  A: {len(ga.g.nodes):5d} nodes  B: {len(gb.g.nodes):5d} nodes")
