#!/usr/bin/env python3
"""
gen_synthetic_gfa.py — Generate synthetic GFA graph pairs for Phase 5 experiments.

Generates two variation graphs (GA, GB) as GFA v1.1 files.
Topology: linear chain of snarls, each snarl a bubble with k alleles.
Allele counts drawn from Zipf distribution (α=2.2) matching HPRC statistics.

Usage:
    # Standard mode (Zipf allele counts)
    python gen_synthetic_gfa.py \\
        --n-snarls 10000 \\
        --out-a ga.gfa --out-b gb.gfa

    # Planted mode (ground truth for Exp1)
    python gen_synthetic_gfa.py \\
        --mode planted \\
        --n-snarls 100 --k-shared 2 --k-only 1 --seed 42 \\
        --out-a ga.gfa --out-b gb.gfa \\
        --truth-a truth_a.txt --truth-b truth_b.txt
"""

import argparse
import random
import sys
import os
from typing import Optional

DNA = "ACGT"


def random_dna(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(DNA) for _ in range(length))


def mutate_seq(seq: str, rate: float, rng: random.Random) -> str:
    return "".join(
        (rng.choice(DNA) if rng.random() < rate else c)
        for c in seq
    )


def zipf_k(rng: random.Random, k_min: int = 2, k_max: int = 50,
           alpha: float = 2.2) -> int:
    """Sample k from truncated Zipf matching HPRC snarl-size distribution."""
    weights = [1.0 / (i ** alpha) for i in range(k_min, k_max + 1)]
    total = sum(weights)
    r = rng.random() * total
    cumulative = 0.0
    for i, w in enumerate(weights):
        cumulative += w
        if r <= cumulative:
            return k_min + i
    return k_max


def write_gfa(nodes: dict, edges: list, paths: dict, out_path: str) -> None:
    """Write GFA v1.1 file."""
    with open(out_path, "w") as f:
        f.write("H\tVN:Z:1.1\n")
        for name, seq in nodes.items():
            f.write(f"S\t{name}\t{seq}\n")
        written_edges = set()
        for u, v in edges:
            key = (u, v)
            if key not in written_edges:
                f.write(f"L\t{u}\t+\t{v}\t+\t0M\n")
                written_edges.add(key)
        for path_name, node_list in paths.items():
            if len(node_list) < 2:
                continue
            segments = ",".join(n + "+" for n in node_list)
            overlaps = ",".join("0M" for _ in range(len(node_list) - 1))
            f.write(f"P\t{path_name}\t{segments}\t{overlaps}\n")


def generate_standard(
    n_snarls: int,
    shared_frac: float = 0.7,
    seq_len_mean: int = 12,
    seq_len_sv: int = 100,
    divergence: float = 0.05,
    seed: int = 42,
    k_max: int = 50,
) -> tuple[dict, dict, dict, dict, list, list]:
    """
    Generate a standard synthetic graph pair with Zipf-distributed allele counts.

    Returns (nodes_a, nodes_b, paths_a, paths_b, edges_a, edges_b)
    """
    rng = random.Random(seed)

    nodes_a, nodes_b = {}, {}
    edges_a, edges_b = [], []
    paths_a_hap, paths_b_hap = [], []

    # Backbone anchors
    anchors = []
    for i in range(n_snarls + 1):
        seq = random_dna(rng.randint(4, 10), rng)
        name = f"anc_{i}"
        nodes_a[name] = seq
        nodes_b[name] = seq
        anchors.append(name)

    # Snarl bubbles
    for s in range(n_snarls):
        k = zipf_k(rng, k_max=k_max)
        n_shared = max(1, round(k * shared_frac))
        n_a_only = rng.randint(0, max(0, k - n_shared))
        n_b_only = max(0, k - n_shared - n_a_only)

        src, snk = anchors[s], anchors[s + 1]
        alleles_a, alleles_b = [], []

        # Shared alleles (same sequence, different node names)
        for j in range(n_shared):
            # Mix of SNP-length and SV-length sequences
            seq_len = seq_len_sv if rng.random() < 0.05 else rng.randint(
                max(1, seq_len_mean // 2), seq_len_mean * 2
            )
            seq = random_dna(seq_len, rng)
            na, nb = f"s{s}_sh{j}_a", f"s{s}_sh{j}_b"
            nodes_a[na] = seq
            nodes_b[nb] = seq
            edges_a += [(src, na), (na, snk)]
            edges_b += [(src, nb), (nb, snk)]
            alleles_a.append(na)
            alleles_b.append(nb)

        # A-only alleles
        for j in range(n_a_only):
            seq_len = seq_len_sv if rng.random() < 0.05 else rng.randint(
                max(1, seq_len_mean // 2), seq_len_mean * 2
            )
            seq = random_dna(seq_len, rng)
            na = f"s{s}_ao{j}"
            nodes_a[na] = seq
            edges_a += [(src, na), (na, snk)]
            alleles_a.append(na)

        # B-only alleles
        for j in range(n_b_only):
            seq_len = seq_len_sv if rng.random() < 0.05 else rng.randint(
                max(1, seq_len_mean // 2), seq_len_mean * 2
            )
            seq = mutate_seq(random_dna(seq_len, rng), divergence, rng)
            nb = f"s{s}_bo{j}"
            nodes_b[nb] = seq
            edges_b += [(src, nb), (nb, snk)]
            alleles_b.append(nb)

        # Ensure minimum 2 alleles per side
        for side, alleles, nodes, pfx in [
            ("a", alleles_a, nodes_a, f"s{s}_fill_a"),
            ("b", alleles_b, nodes_b, f"s{s}_fill_b"),
        ]:
            while len(alleles) < 2:
                seq = random_dna(seq_len_mean, rng)
                name = f"{pfx}{len(alleles)}"
                nodes[name] = seq
                if side == "a":
                    edges_a += [(src, name), (name, snk)]
                else:
                    edges_b += [(src, name), (name, snk)]
                alleles.append(name)

        # Representative haplotype paths (allele[0])
        paths_a_hap.append(alleles_a[0])
        paths_b_hap.append(alleles_b[0])

    # Stitch haplotype paths
    path_a = [anchors[0]]
    path_b = [anchors[0]]
    for s in range(n_snarls):
        path_a += [paths_a_hap[s], anchors[s + 1]]
        path_b += [paths_b_hap[s], anchors[s + 1]]

    paths_a = {"hap_A0": path_a}
    paths_b = {"hap_B0": path_b}

    return nodes_a, nodes_b, paths_a, paths_b, edges_a, edges_b


def generate_planted(
    n_snarls: int,
    k_shared: int = 2,
    k_only: int = 1,
    seq_len: int = 12,
    seed: int = 42,
) -> tuple[dict, dict, dict, dict, list, list, list, list]:
    """
    Generate a planted-variant graph pair for Exp1 correctness evaluation.

    Returns (nodes_a, nodes_b, paths_a, paths_b, edges_a, edges_b,
             truth_a_seqs, truth_b_seqs)
    """
    rng = random.Random(seed)

    nodes_a, nodes_b = {}, {}
    edges_a, edges_b = [], []
    truth_a, truth_b = [], []

    anchors = []
    for i in range(n_snarls + 1):
        seq = random_dna(8, rng)
        name = f"anc_{i}"
        nodes_a[name] = seq
        nodes_b[name] = seq
        anchors.append(name)

    paths_a_hap, paths_b_hap = [], []

    for s in range(n_snarls):
        src, snk = anchors[s], anchors[s + 1]
        alleles_a, alleles_b = [], []

        for j in range(k_shared):
            seq = random_dna(seq_len, rng)
            na, nb = f"s{s}_sh{j}_a", f"s{s}_sh{j}_b"
            nodes_a[na] = seq
            nodes_b[nb] = seq
            edges_a += [(src, na), (na, snk)]
            edges_b += [(src, nb), (nb, snk)]
            alleles_a.append(na)
            alleles_b.append(nb)

        for j in range(k_only):
            seq = random_dna(seq_len, rng)
            truth_a.append(seq)
            na = f"s{s}_ao{j}"
            nodes_a[na] = seq
            edges_a += [(src, na), (na, snk)]
            alleles_a.append(na)

        for j in range(k_only):
            seq = random_dna(seq_len, rng)
            truth_b.append(seq)
            nb = f"s{s}_bo{j}"
            nodes_b[nb] = seq
            edges_b += [(src, nb), (nb, snk)]
            alleles_b.append(nb)

        paths_a_hap.append(alleles_a[0])
        paths_b_hap.append(alleles_b[0])

    path_a = [anchors[0]]
    path_b = [anchors[0]]
    for s in range(n_snarls):
        path_a += [paths_a_hap[s], anchors[s + 1]]
        path_b += [paths_b_hap[s], anchors[s + 1]]

    return (nodes_a, nodes_b,
            {"hap_A0": path_a}, {"hap_B0": path_b},
            edges_a, edges_b,
            truth_a, truth_b)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", choices=["standard", "planted"], default="standard")
    parser.add_argument("--n-snarls", type=int, default=1000)
    parser.add_argument("--shared-frac", type=float, default=0.7)
    parser.add_argument("--seq-len", type=int, default=12)
    parser.add_argument("--seq-len-sv", type=int, default=100,
                        help="Sequence length for SV alleles (5%% of snarls)")
    parser.add_argument("--divergence", type=float, default=0.05)
    parser.add_argument("--k-max", type=int, default=50,
                        help="Max alleles per snarl (Zipf truncation)")
    parser.add_argument("--k-shared", type=int, default=2,
                        help="Shared alleles per snarl (planted mode)")
    parser.add_argument("--k-only", type=int, default=1,
                        help="A-only / B-only alleles per snarl (planted mode)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-a", required=True)
    parser.add_argument("--out-b", required=True)
    parser.add_argument("--truth-a", default=None,
                        help="Write planted A-only sequences (one per line)")
    parser.add_argument("--truth-b", default=None,
                        help="Write planted B-only sequences (one per line)")
    args = parser.parse_args()

    if args.mode == "planted":
        na, nb, pa, pb, ea, eb, ta, tb = generate_planted(
            n_snarls=args.n_snarls,
            k_shared=args.k_shared,
            k_only=args.k_only,
            seq_len=args.seq_len,
            seed=args.seed,
        )
        if args.truth_a:
            with open(args.truth_a, "w") as f:
                f.write("\n".join(ta) + "\n")
        if args.truth_b:
            with open(args.truth_b, "w") as f:
                f.write("\n".join(tb) + "\n")
    else:
        na, nb, pa, pb, ea, eb = generate_standard(
            n_snarls=args.n_snarls,
            shared_frac=args.shared_frac,
            seq_len_mean=args.seq_len,
            seq_len_sv=args.seq_len_sv,
            divergence=args.divergence,
            seed=args.seed,
            k_max=args.k_max,
        )

    write_gfa(na, ea, pa, args.out_a)
    write_gfa(nb, eb, pb, args.out_b)

    a_nodes = len(na)
    b_nodes = len(nb)
    a_size = os.path.getsize(args.out_a)
    b_size = os.path.getsize(args.out_b)

    print(f"[gen] mode={args.mode}  n_snarls={args.n_snarls}  seed={args.seed}",
          file=sys.stderr)
    print(f"[gen] GA: {a_nodes} nodes  {a_size/1024:.1f} KB → {args.out_a}",
          file=sys.stderr)
    print(f"[gen] GB: {b_nodes} nodes  {b_size/1024:.1f} KB → {args.out_b}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
