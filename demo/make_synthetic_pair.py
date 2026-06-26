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
import random

DNA = "ACGT"


def write_gfa(path, ref, variants, chop, label, ref_name):
    """variants: {pos0 -> alt_base}. Reference is chopped into `chop`-sized
    nodes, always breaking at variant positions so every bubble opens at the
    variant's exact reference coordinate regardless of `chop`."""
    nodes, edges = [], []
    nid = [0]

    def add(seq):
        name = f"{label}{nid[0]}"
        nid[0] += 1
        nodes.append((name, seq))
        return name

    backbone = []  # reference-path nodes: {name, is_var, alt}
    i = 0
    for p in sorted(variants):
        run = ref[i:p]                                  # invariant stretch
        for c in range(0, len(run), chop):
            backbone.append({"name": add(run[c:c + chop]), "is_var": False, "alt": None})
        backbone.append({"name": add(ref[p]),           # 1bp reference allele
                         "is_var": True,
                         "alt": add(variants[p])})       # 1bp alternate allele
        i = p + 1
    run = ref[i:]
    for c in range(0, len(run), chop):
        backbone.append({"name": add(run[c:c + chop]), "is_var": False, "alt": None})

    for k in range(len(backbone) - 1):                  # reference-path edges
        edges.append((backbone[k]["name"], backbone[k + 1]["name"]))
    for k, b in enumerate(backbone):                    # alt bubble edges
        if b["is_var"]:
            assert 0 < k < len(backbone) - 1, "variant too close to an end"
            edges.append((backbone[k - 1]["name"], b["alt"]))
            edges.append((b["alt"], backbone[k + 1]["name"]))

    with open(path, "w") as f:
        f.write("H\tVN:Z:1.1\n")
        for name, seq in nodes:
            f.write(f"S\t{name}\t{seq}\n")
        seen = set()
        for u, v in edges:
            if (u, v) not in seen:
                f.write(f"L\t{u}\t+\t{v}\t+\t0M\n")
                seen.add((u, v))
        segs = ",".join(n["name"] + "+" for n in backbone)
        f.write(f"P\t{ref_name}\t{segs}\t*\n")           # reference path (P-line)


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
    a = ap.parse_args()

    rng = random.Random(a.seed)
    ref = "".join(rng.choice(DNA) for _ in range(a.length))

    margin = max(a.chop_a, a.chop_b) + 2
    pool = list(range(margin, a.length - margin))
    rng.shuffle(pool)
    chosen = []
    for c in sorted(pool):                 # keep variants >=2 apart (clean bubbles)
        if all(abs(c - x) >= 2 for x in chosen):
            chosen.append(c)
    need = a.n_shared + a.n_a + a.n_b
    assert len(chosen) >= need, "increase --length or lower variant counts"
    chosen = chosen[:need]
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

    varsA = {p: altof[p] for p in cls if cls[p] in ("shared", "A")}
    varsB = {p: altof[p] for p in cls if cls[p] in ("shared", "B")}

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
