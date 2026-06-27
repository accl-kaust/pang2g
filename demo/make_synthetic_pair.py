#!/usr/bin/env python3
"""
make_synthetic_pair.py — controlled reference-anchored pangenome pair.

Builds two GFA graphs (A and B) that embed the SAME reference path but plant a
known set of variants and — crucially — chop the shared reference into
DIFFERENT nodes.  This is the situation that breaks naive boundary matching and
that `g2g --ref` is designed to handle: variants are matched by reference
COORDINATE, not by node boundaries.

Plants:
  * n_shared SNPs present (same alt) in BOTH graphs  -> should be "shared"
  * n_a      SNPs present only in A                   -> should be A-only
  * n_b      SNPs present only in B                   -> should be B-only

Outputs demoA.gfa, demoB.gfa and truth.tsv (the ground truth).
"""
import argparse
import os
import random

DNA = "ACGT"


def write_gfa(path, ref, variants, chop, label, ref_name):
    """Stream a GFA for `ref` with a 1 bp bubble planted at each `variants`
    position ({pos0 -> alt_base}). The reference is chopped into `chop`-sized
    nodes, always breaking at variant positions so every bubble opens at the
    variant's exact reference coordinate regardless of `chop`.

    Written in O(1) memory via two deterministic walks (the P-line, then the
    S/L lines) instead of materialising every node/edge — so the generator
    scales to hundreds of megabases. g2g's GFA reader is two-pass / order
    independent, so emitting P before S/L is fine.
    """
    positions = sorted(variants)

    def walk():
        """Yield (kind, seq) for every node in canonical creation order;
        kind in {'run','ref','alt'}. The name counter advances once per yield,
        so both passes below produce identical node names."""
        i = 0
        for p in positions:
            run = ref[i:p]                              # invariant stretch
            for c in range(0, len(run), chop):
                yield "run", run[c:c + chop]
            yield "ref", ref[p]                         # reference allele
            alts = variants[p]                          # one alt (str) or a list
            for alt in (alts,) if isinstance(alts, str) else alts:
                yield "alt", alt                        # 1+ alternate alleles
            i = p + 1
        run = ref[i:]
        for c in range(0, len(run), chop):
            yield "run", run[c:c + chop]

    with open(path, "w") as f:
        f.write("H\tVN:Z:1.1\n")

        # Pass 1 — reference P-line (backbone nodes only, in order).
        f.write(f"P\t{ref_name}\t")
        nid = 0
        first = True
        for kind, _ in walk():
            nm = f"{label}{nid}"; nid += 1
            if kind != "alt":
                f.write(f"{nm}+" if first else f",{nm}+")
                first = False
        f.write("\t*\n")

        # Pass 2 — S (segments) and L (links), streamed.
        nid = 0
        prev_bb = None        # previous backbone node (reference-path chain)
        b_before = None       # backbone node just before the current variant
        pending = []          # alt nodes awaiting their closing edge
        for kind, seq in walk():
            nm = f"{label}{nid}"; nid += 1
            f.write(f"S\t{nm}\t{seq}\n")
            if kind == "alt":
                if b_before is not None:               # bubble opens
                    f.write(f"L\t{b_before}\t+\t{nm}\t+\t0M\n")
                pending.append(nm)
                continue
            # backbone node ('run' or 'ref'): extend chain, close any open bubble
            if prev_bb is not None:
                f.write(f"L\t{prev_bb}\t+\t{nm}\t+\t0M\n")
            for a in pending:                          # bubble closes here
                f.write(f"L\t{a}\t+\t{nm}\t+\t0M\n")
            pending = []
            if kind == "ref":
                b_before = prev_bb
            prev_bb = nm


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--length", type=int, default=400)
    ap.add_argument("--n-shared", type=int, default=5)
    ap.add_argument("--n-a", type=int, default=3)
    ap.add_argument("--n-b", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--chop-a", type=int, default=25, help="node size in graph A")
    ap.add_argument("--chop-b", type=int, default=17, help="node size in graph B")
    ap.add_argument("--ref-name", default="GRCh38")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--max-k", type=int, default=2,
                    help="max alleles per snarl; >2 plants heavy-tailed "
                         "multi-allele snarls (work for the parallel aligner)")
    ap.add_argument("--complex-frac", type=float, default=0.0,
                    help="fraction of shared sites made multi-allele when --max-k>2")
    ap.add_argument("--alt-len", type=int, default=4,
                    help="max length of multi-allele alt sequences (long = SV-like, "
                         "heavy O(L^2) per-pair edit distance for the aligner DP)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    rng = random.Random(a.seed)
    # Random reference, generated in chunks so multi-megabase references stay
    # fast and memory-light (a per-base Python loop dominated runtime, and a
    # single rng.choices(k=length) would build a huge intermediate list).
    ref_ba = bytearray(a.length)
    pos, CHUNK = 0, 1 << 20
    while pos < a.length:
        k = min(CHUNK, a.length - pos)
        ref_ba[pos:pos + k] = bytes(rng.choices(b"ACGT", k=k))
        pos += k
    ref = ref_ba.decode("ascii")
    del ref_ba

    margin = max(a.chop_a, a.chop_b) + 2
    need = a.n_shared + a.n_a + a.n_b
    usable = a.length - 2 * margin
    # Place variants on an even grid so they are always >=2 bp apart (clean,
    # independent bubbles) regardless of count, and spread across the whole
    # reference. O(N) — scales to multi-megabase references (the earlier
    # shuffle-pool + pairwise proximity check was O(length * need)).
    assert usable >= need * 2, "increase --length or lower variant counts"
    step = usable // need
    chosen = [margin + i * step for i in range(need)]
    rng.shuffle(chosen)
    shared = chosen[:a.n_shared]
    aonly = chosen[a.n_shared:a.n_shared + a.n_a]
    bonly = chosen[a.n_shared + a.n_a:]

    def alt(p):
        return rng.choice([b for b in DNA if b != ref[p]])

    cls = {}
    for p in shared: cls[p] = "shared"
    for p in aonly:  cls[p] = "A"
    for p in bonly:  cls[p] = "B"
    altof = {p: alt(p) for p in cls}

    # Optionally upgrade a fraction of SHARED sites into heavy-tailed multi-allele
    # snarls (k alleles, k drawn Pareto and capped at --max-k). These give the
    # work-stealing aligner real, skewed per-snarl matching work; the --ref diff
    # is unaffected because A and B carry the same allele set at these sites.
    complex_alts = {}
    if a.max_k > 2 and a.complex_frac > 0:
        for p in shared:
            if rng.random() < a.complex_frac:
                k = min(a.max_k, 2 + int(rng.paretovariate(1.3)))
                alts = set()
                while len(alts) < k - 1:
                    L = rng.randint(1, a.alt_len)
                    alts.add("".join(rng.choice(DNA) for _ in range(L)))
                complex_alts[p] = sorted(alts)

    def alleles(p):
        return complex_alts[p] if p in complex_alts else altof[p]

    varsA = {p: alleles(p) for p in cls if cls[p] in ("shared", "A")}
    varsB = {p: alleles(p) for p in cls if cls[p] in ("shared", "B")}

    write_gfa(f"{a.outdir}/demoA.gfa", ref, varsA, a.chop_a, "A", a.ref_name)
    write_gfa(f"{a.outdir}/demoB.gfa", ref, varsB, a.chop_b, "B", a.ref_name)
    with open(f"{a.outdir}/truth.tsv", "w") as f:
        f.write("pos\tref\talt\tclass\n")
        for p in sorted(cls):
            f.write(f"{p}\t{ref[p]}\t{altof[p]}\t{cls[p]}\n")

    print(f"[gen] ref={a.length}bp  shared={a.n_shared} A_only={a.n_a} "
          f"B_only={a.n_b}  chop A={a.chop_a} B={a.chop_b} (different on purpose)")


if __name__ == "__main__":
    main()
