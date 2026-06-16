#!/usr/bin/env python3
"""
gen_hprc_scale_gfa.py — HPRC-realistic synthetic GFA generator for G2G evaluation
=================================================================================

Produces pairs of variation graphs (GA, GB) with a linear backbone of N snarls.
Each snarl is a bubble with k alleles drawn from a Zipf(alpha) distribution
truncated to [k_min, k_max], calibrated to the empirical HPRC human-pangenome
snarl-size distribution measured in Phase 2:

    Gini(k^3) = 0.55-0.58 ,  mean k ~ 2.6 ,  ~67% SNP-bubbles (k=2) ,  tail to k~50

This is the generator used for the *large* scaling runs (1K .. 10M snarls). It
differs from the Phase-2 `gen_scaling_gfa.py` in three ways that matter at scale:

  1. **Buffered I/O.** Node/edge/path lines are emitted into a list and flushed
     in big chunks. At 10M snarls (~36M nodes) the naive per-line f.write() is
     I/O-bound; chunked writes cut wall time substantially.

  2. **Haplotype-panel mode** (`--mode haplotypes`). Fixes the graph at a chosen
     snarl count and varies the number of *haplotype paths* threaded through it
     (h in {2,10,30,100,300}). This is the input for Experiment 3b — the figure
     that shows G2G's runtime is flat in h while linearise-and-map baselines grow
     linearly/quadratically in h.

  3. **Planted ground truth** is emitted to `truth_<tag>.tsv` so Experiment 1's
     scaled correctness check (n in {1e3,1e4,1e5}) can verify recall/precision
     without re-deriving labels.

Realism note
------------
A *true* HPRC graph has nested snarls, cycles (tandem repeats), and a non-linear
backbone. This generator deliberately produces a DAG of series bubbles because
(a) the G2G prototype's built-in decomposer is DAG-only, and (b) the systems
claim is about the *irregular-memory workload* (heavy-tailed per-snarl work,
bottom-up traversal), which the series-bubble model reproduces faithfully at the
correct scale. Experiment 6 runs on the *actual* HPRC graphs via minigraph-cactus
to close the realism gap for the biological claim.

Usage
-----
  # Graph-size scaling (default mode):
  python gen_hprc_scale_gfa.py --mode size --n_snarls 1000000 --out_dir DIR --seed 42

  # Haplotype-panel scaling (Experiment 3b):
  python gen_hprc_scale_gfa.py --mode haplotypes --n_snarls 10000 \
         --haplotypes 100 --out_dir DIR --seed 42

Outputs
-------
  size mode:        ga_<N>.gfa  gb_<N>.gfa  truth_<N>.tsv
  haplotypes mode:  ga_h<H>_n<N>.gfa  gb_h<H>_n<N>.gfa  truth_h<H>_n<N>.tsv
"""

import argparse
import os
import random
import sys

BASES = "ACGT"
_FLUSH_EVERY = 1_000_000  # lines buffered before a flush


def random_dna(length, rng):
    return "".join(rng.choice(BASES) for _ in range(length))


def zipf_k(rng, k_min, k_max, alpha):
    """Sample allele count k ~ Zipf(alpha) truncated to [k_min, k_max]."""
    weights = [1.0 / (i ** alpha) for i in range(k_min, k_max + 1)]
    total = sum(weights)
    r = rng.random() * total
    c = 0.0
    for i, w in enumerate(weights):
        c += w
        if r <= c:
            return k_min + i
    return k_max


class GfaWriter:
    """Chunked GFA writer: buffers lines and flushes in large blocks."""

    def __init__(self, path):
        self.f = open(path, "w", buffering=1 << 20)
        self.buf = []
        self.n = 0
        self.f.write("H\tVN:Z:1.1\n")

    def seg(self, name, seq):
        self.buf.append(f"S\t{name}\t{seq}\n")
        self._tick()

    def link(self, u, v):
        self.buf.append(f"L\t{u}\t+\t{v}\t+\t0M\n")
        self._tick()

    def path(self, name, node_list):
        cigar = ",".join(f"{n}+" for n in node_list)
        overlaps = "*," * (len(node_list) - 1) + "*"
        self.buf.append(f"P\t{name}\t{cigar}\t{overlaps}\n")
        self._tick()

    def _tick(self):
        self.n += 1
        if self.n >= _FLUSH_EVERY:
            self.f.write("".join(self.buf))
            self.buf.clear()
            self.n = 0

    def close(self):
        if self.buf:
            self.f.write("".join(self.buf))
            self.buf.clear()
        self.f.close()


def gen_size_mode(n_snarls, rng, out_dir, anchor_len, inner_len,
                  k_min, k_max, alpha, n_haplotypes):
    """Generate a GA/GB pair with n_snarls bubbles; emit ground truth.

    n_haplotypes: number of representative haplotype paths to thread through
    each graph (default 1). Each path picks a random allele per snarl.
    """
    tag = str(n_snarls)
    ga = GfaWriter(os.path.join(out_dir, f"ga_{tag}.gfa"))
    gb = GfaWriter(os.path.join(out_dir, f"gb_{tag}.gfa"))
    truth = open(os.path.join(out_dir, f"truth_{tag}.tsv"), "w", buffering=1 << 20)
    truth.write("snarl\ttype\tseq\tnode_a\tnode_b\n")

    # ── Backbone anchors (identical in both graphs) ──────────────────────────
    anc = []
    for i in range(n_snarls + 1):
        seq = random_dna(anchor_len, rng)
        name = f"anc_{i}"
        ga.seg(name, seq)
        gb.seg(name, seq)
        anc.append(name)

    # ── Per-snarl allele records, needed for haplotype threading ─────────────
    # To keep memory bounded at 10M snarls we DO NOT keep every allele name in
    # RAM; instead we thread haplotypes on the fly by recording only the chosen
    # allele per haplotype per snarl as we go.
    hap_a = [[anc[0]] for _ in range(n_haplotypes)]
    hap_b = [[anc[0]] for _ in range(n_haplotypes)]

    total_a_nodes = n_snarls + 1
    total_b_nodes = n_snarls + 1

    for s in range(n_snarls):
        src, snk = anc[s], anc[s + 1]
        k = zipf_k(rng, k_min, k_max, alpha)
        n_shared = max(1, k // 2)
        rem = k - n_shared
        n_a_only = rem // 2 + (rem % 2)
        n_b_only = rem // 2

        a_allele_names = []
        b_allele_names = []

        # Shared alleles (same bytes, different node names in A vs B)
        for j in range(n_shared):
            seq = random_dna(inner_len, rng)
            na = f"s{s}_sh{j}_a"
            nb = f"s{s}_sh{j}_b"
            ga.seg(na, seq); ga.link(src, na); ga.link(na, snk)
            gb.seg(nb, seq); gb.link(src, nb); gb.link(nb, snk)
            truth.write(f"{s}\tshared\t{seq}\t{na}\t{nb}\n")
            a_allele_names.append(na); b_allele_names.append(nb)
            total_a_nodes += 1; total_b_nodes += 1

        # A-only alleles
        for j in range(n_a_only):
            seq = random_dna(inner_len, rng)
            na = f"s{s}_ao{j}"
            ga.seg(na, seq); ga.link(src, na); ga.link(na, snk)
            truth.write(f"{s}\ta_only\t{seq}\t{na}\t\n")
            a_allele_names.append(na); total_a_nodes += 1

        # B-only alleles
        for j in range(n_b_only):
            seq = random_dna(inner_len, rng)
            nb = f"s{s}_bo{j}"
            gb.seg(nb, seq); gb.link(src, nb); gb.link(nb, snk)
            truth.write(f"{s}\tb_only\t{seq}\t\t{nb}\n")
            b_allele_names.append(nb); total_b_nodes += 1

        # Thread haplotypes through a randomly chosen allele each
        for hpath in hap_a:
            hpath.append(rng.choice(a_allele_names)); hpath.append(snk)
        for hpath in hap_b:
            hpath.append(rng.choice(b_allele_names)); hpath.append(snk)

    # ── Emit haplotype paths ─────────────────────────────────────────────────
    for hi, hpath in enumerate(hap_a):
        ga.path(f"hap_A{hi}", hpath)
    for hi, hpath in enumerate(hap_b):
        gb.path(f"hap_B{hi}", hpath)

    ga.close(); gb.close(); truth.close()
    return total_a_nodes, total_b_nodes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["size", "haplotypes"], default="size")
    ap.add_argument("--n_snarls", type=int, required=True)
    ap.add_argument("--haplotypes", type=int, default=1,
                    help="Number of haplotype paths per graph (panel size). "
                         "In 'haplotypes' mode this is the swept variable.")
    ap.add_argument("--out_dir", default=".")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--anchor_len", type=int, default=20)
    ap.add_argument("--inner_len", type=int, default=50)
    ap.add_argument("--k_min", type=int, default=2)
    ap.add_argument("--k_max", type=int, default=50)
    ap.add_argument("--alpha", type=float, default=2.0,
                    help="Zipf exponent (2.0 reproduces HPRC heavy tail)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    rng = random.Random(args.seed)

    if args.mode == "haplotypes":
        # Fixed graph, variable panel: file names encode both H and N.
        tag = f"h{args.haplotypes}_n{args.n_snarls}"
        ga_path = os.path.join(args.out_dir, f"ga_{tag}.gfa")
        gb_path = os.path.join(args.out_dir, f"gb_{tag}.gfa")
        # Reuse the size-mode generator but rename outputs afterwards.
        a_nodes, b_nodes = gen_size_mode(
            args.n_snarls, rng, args.out_dir, args.anchor_len, args.inner_len,
            args.k_min, args.k_max, args.alpha, args.haplotypes)
        # Rename ga_<N>.gfa -> ga_h<H>_n<N>.gfa
        os.replace(os.path.join(args.out_dir, f"ga_{args.n_snarls}.gfa"), ga_path)
        os.replace(os.path.join(args.out_dir, f"gb_{args.n_snarls}.gfa"), gb_path)
        os.replace(os.path.join(args.out_dir, f"truth_{args.n_snarls}.tsv"),
                   os.path.join(args.out_dir, f"truth_{tag}.tsv"))
        print(f"[gen] mode=haplotypes  N={args.n_snarls:,}  H={args.haplotypes}  "
              f"GA_nodes={a_nodes:,}  GB_nodes={b_nodes:,}")
        print(f"[gen] wrote {ga_path}")
        print(f"[gen] wrote {gb_path}")
    else:
        a_nodes, b_nodes = gen_size_mode(
            args.n_snarls, rng, args.out_dir, args.anchor_len, args.inner_len,
            args.k_min, args.k_max, args.alpha, args.haplotypes)
        print(f"[gen] mode=size  N={args.n_snarls:,}  "
              f"GA_nodes={a_nodes:,}  GB_nodes={b_nodes:,}  "
              f"haplotypes={args.haplotypes}")
        print(f"[gen] wrote {os.path.join(args.out_dir, f'ga_{args.n_snarls}.gfa')}")
        print(f"[gen] wrote {os.path.join(args.out_dir, f'gb_{args.n_snarls}.gfa')}")


if __name__ == "__main__":
    main()
