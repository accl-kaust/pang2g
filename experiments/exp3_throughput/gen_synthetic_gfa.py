#!/usr/bin/env python3
"""
gen_synthetic_gfa.py — Generate synthetic GFA graph pairs for the experiments.

Generates two variation graphs (GA, GB) as GFA v1.1 files.
Topology: linear chain of snarls, each snarl a bubble with k alleles.
Allele counts drawn from Zipf distribution (alpha=2.2) matching HPRC statistics.

Usage:
    # Standard mode (Zipf allele counts; throughput / scaling / ablation inputs)
    python gen_synthetic_gfa.py \\
        --n-snarls 10000 \\
        --out-a ga.gfa --out-b gb.gfa

    # Planted mode (ground truth for Exp1: per-snarl shared / A-only / B-only alleles)
    python gen_synthetic_gfa.py \\
        --mode planted \\
        --n-snarls 100 --k-shared 2 --k-only 1 --seed 42 \\
        --out-a ga.gfa --out-b gb.gfa \\
        --truth-a truth_a.txt --truth-b truth_b.txt

Planted-mode knobs that stress snarl PAIRING (Exp1b, pairing error vs. size):
    --anchor-len L      length of the snarl boundary nodes (default 8).  Short
                        boundaries (1-2 bp) make boundary-sequence keys collide,
                        which is what real vg-construct graphs look like next to
                        dense variation.
    --private-frac f    fraction of snarl slots that exist in only ONE graph
                        (alternating A-private / B-private).  The other graph
                        has a plain backbone edge there.  A private snarl shifts
                        the positional (FIFO) pairing of every later snarl that
                        shares its boundary key.  All alleles of a private snarl
                        are ground-truth delta alleles of the graph that has it.
    --linker-every m    insert an anchor-able linker node (in/out-degree 1,
                        --linker-len bp, identical in A and B) after every m-th
                        snarl slot, so the two graphs share anchors between
                        groups of m snarls.  0 = no linkers (default).
    --linker-len L      linker node length (default 24; >= the minimizer k=15).
"""

import argparse
import random
import sys
import os

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
) -> tuple:
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
    anchor_len: int = 8,
    private_frac: float = 0.0,
    linker_every: int = 0,
    linker_len: int = 24,
) -> tuple:
    """
    Generate a planted-variant graph pair for Exp1 correctness evaluation.

    Every snarl slot s has a source boundary node and a sink boundary node of
    `anchor_len` bp (identical sequence in A and B).  A *shared* slot carries
    k_shared alleles present in both graphs plus k_only A-only and k_only
    B-only alleles.  A *private* slot (chosen with probability `private_frac`,
    alternating A/B) carries k_shared + k_only alleles in one graph and a plain
    backbone edge in the other; all of its alleles are delta alleles of the
    graph that has it.

    Returns (nodes_a, nodes_b, paths_a, paths_b, edges_a, edges_b,
             truth_a_seqs, truth_b_seqs, stats)
    """
    rng = random.Random(seed)

    nodes_a, nodes_b = {}, {}
    edges_a, edges_b = [], []
    truth_a, truth_b = [], []
    stats = {"shared_slots": 0, "a_private_slots": 0, "b_private_slots": 0,
             "linkers": 0}

    # Slot kinds: 'S' shared, 'A' A-private, 'B' B-private.
    kinds = []
    flip = 0
    for s in range(n_snarls):
        if private_frac > 0 and rng.random() < private_frac:
            kinds.append("A" if flip == 0 else "B")
            flip ^= 1
        else:
            kinds.append("S")

    # Boundary nodes.  Without a linker after slot s, the sink of slot s IS the
    # source of slot s+1 (a shared boundary, as in a dense variant run); with a
    # linker, sink -> lnk_s -> next source, and the linker is an anchor.
    path_a, path_b = [], []
    src_name = "b_0"
    nodes_a[src_name] = nodes_b[src_name] = random_dna(anchor_len, rng)
    path_a.append(src_name)
    path_b.append(src_name)

    for s in range(n_snarls):
        src = src_name
        has_linker = (linker_every > 0 and (s + 1) % linker_every == 0
                      and s + 1 < n_snarls)
        snk = f"e_{s}" if has_linker else f"b_{s + 1}"
        nodes_a[snk] = nodes_b[snk] = random_dna(anchor_len, rng)

        kind = kinds[s]
        alleles_a, alleles_b = [], []

        if kind == "S":
            stats["shared_slots"] += 1
            for j in range(k_shared):
                aseq = random_dna(seq_len, rng)
                na, nb = f"s{s}_sh{j}_a", f"s{s}_sh{j}_b"
                nodes_a[na] = aseq
                nodes_b[nb] = aseq
                edges_a += [(src, na), (na, snk)]
                edges_b += [(src, nb), (nb, snk)]
                alleles_a.append(na)
                alleles_b.append(nb)
            for j in range(k_only):
                aseq = random_dna(seq_len, rng)
                truth_a.append(aseq)
                na = f"s{s}_ao{j}"
                nodes_a[na] = aseq
                edges_a += [(src, na), (na, snk)]
                alleles_a.append(na)
            for j in range(k_only):
                bseq = random_dna(seq_len, rng)
                truth_b.append(bseq)
                nb = f"s{s}_bo{j}"
                nodes_b[nb] = bseq
                edges_b += [(src, nb), (nb, snk)]
                alleles_b.append(nb)
        else:
            # Private slot: bubble in one graph only; the other graph walks the
            # backbone edge src -> snk directly.
            if kind == "A":
                own_nodes, own_edges, own_truth, own_alleles, tag = (
                    nodes_a, edges_a, truth_a, alleles_a, "ap")
                stats["a_private_slots"] += 1
                edges_b.append((src, snk))
            else:
                own_nodes, own_edges, own_truth, own_alleles, tag = (
                    nodes_b, edges_b, truth_b, alleles_b, "bp")
                stats["b_private_slots"] += 1
                edges_a.append((src, snk))
            for j in range(k_shared + k_only):
                pseq = random_dna(seq_len, rng)
                own_truth.append(pseq)
                nm = f"s{s}_{tag}{j}"
                own_nodes[nm] = pseq
                own_edges += [(src, nm), (nm, snk)]
                own_alleles.append(nm)

        # Haplotype paths (first allele where a bubble exists, else direct)
        if alleles_a:
            path_a.append(alleles_a[0])
        path_a.append(snk)
        if alleles_b:
            path_b.append(alleles_b[0])
        path_b.append(snk)

        if has_linker:
            lnk = f"lnk_{s}"
            nodes_a[lnk] = nodes_b[lnk] = random_dna(linker_len, rng)
            nxt = f"b_{s + 1}"
            nodes_a[nxt] = nodes_b[nxt] = random_dna(anchor_len, rng)
            edges_a += [(snk, lnk), (lnk, nxt)]
            edges_b += [(snk, lnk), (lnk, nxt)]
            path_a += [lnk, nxt]
            path_b += [lnk, nxt]
            stats["linkers"] += 1
            src_name = nxt
        else:
            src_name = snk

    return (nodes_a, nodes_b,
            {"hap_A0": path_a}, {"hap_B0": path_b},
            edges_a, edges_b,
            truth_a, truth_b, stats)


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
    parser.add_argument("--anchor-len", type=int, default=8,
                        help="Boundary node length in bp (planted mode)")
    parser.add_argument("--private-frac", type=float, default=0.0,
                        help="Fraction of snarl slots private to one graph (planted mode)")
    parser.add_argument("--linker-every", type=int, default=0,
                        help="Insert an anchor-able linker node after every m-th slot (planted mode)")
    parser.add_argument("--linker-len", type=int, default=24,
                        help="Linker node length in bp (planted mode)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-a", required=True)
    parser.add_argument("--out-b", required=True)
    parser.add_argument("--truth-a", default=None,
                        help="Write planted A-only sequences (one per line)")
    parser.add_argument("--truth-b", default=None,
                        help="Write planted B-only sequences (one per line)")
    args = parser.parse_args()

    stats = {}
    if args.mode == "planted":
        na, nb, pa, pb, ea, eb, ta, tb, stats = generate_planted(
            n_snarls=args.n_snarls,
            k_shared=args.k_shared,
            k_only=args.k_only,
            seq_len=args.seq_len,
            seed=args.seed,
            anchor_len=args.anchor_len,
            private_frac=args.private_frac,
            linker_every=args.linker_every,
            linker_len=args.linker_len,
        )
        if args.truth_a:
            with open(args.truth_a, "w") as f:
                f.write("\n".join(ta) + ("\n" if ta else ""))
        if args.truth_b:
            with open(args.truth_b, "w") as f:
                f.write("\n".join(tb) + ("\n" if tb else ""))
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

    print(f"[gen] mode={args.mode}  n_snarls={args.n_snarls}  seed={args.seed}"
          + (f"  slots={stats}" if stats else ""),
          file=sys.stderr)
    print(f"[gen] GA: {a_nodes} nodes  {a_size/1024:.1f} KB -> {args.out_a}",
          file=sys.stderr)
    print(f"[gen] GB: {b_nodes} nodes  {b_size/1024:.1f} KB -> {args.out_b}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
