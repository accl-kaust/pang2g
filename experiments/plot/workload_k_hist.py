#!/usr/bin/env python3
"""
workload_k_hist.py - allele-count (k) histogram of the snarls in a GFA pair.

Reads two GFA files written by the project's generators (simple bubbles:
source -> {alleles} -> sink), reports the k distribution of
  * all snarls of A, all snarls of B, and
  * the snarls that the boundary-key pairing (correspond_snarls) pairs,
and writes a CSV (k,count) for make_figures.py --alleles-csv.

Usage:
    python3 experiments/plot/workload_k_hist.py A.gfa B.gfa --out k_hist.csv
"""
import argparse
import csv
from collections import Counter, defaultdict, deque


def read_gfa(path):
    seq, succ, pred = {}, defaultdict(list), defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("S\t"):
                p = line.rstrip("\n").split("\t")
                seq[p[1]] = p[2]
            elif line.startswith("L\t"):
                p = line.rstrip("\n").split("\t")
                succ[p[1]].append(p[3])
                pred[p[3]].append(p[1])
    return seq, succ, pred


def bubbles(seq, succ):
    """[(src, sink, k)] for simple bubbles src -> alleles -> sink, in node order
    of the source names (the generators name nodes in backbone order)."""
    out = []
    for src, outs in succ.items():
        if len(outs) < 2:
            continue
        sinks = Counter()
        for a in outs:
            nxt = succ.get(a, [])
            if len(nxt) == 1:
                sinks[nxt[0]] += 1
            elif a in succ and not nxt:
                pass
        # direct edge src->sink (deletion / private-slot backbone edge)
        for a in outs:
            if a in sinks:
                sinks[a] += 1
        if not sinks:
            continue
        sink, n = sinks.most_common(1)[0]
        out.append((src, sink, len(outs)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gfa_a"); ap.add_argument("gfa_b")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    sa, ta, _ = read_gfa(a.gfa_a)
    sb, tb, _ = read_gfa(a.gfa_b)
    ba, bb = bubbles(sa, ta), bubbles(sb, tb)

    def key(seq, src, sink):
        return seq[src] + "\x01" + seq[sink]

    # replicate correspond_snarls: FIFO within identical boundary keys
    byk = defaultdict(deque)
    for src, sink, k in bb:
        byk[key(sb, src, sink)].append(k)
    paired = []
    for src, sink, k in ba:
        q = byk.get(key(sa, src, sink))
        if q:
            paired.append((k, q.popleft()))

    def report(label, ks):
        c = Counter(ks); n = sum(c.values())
        k2 = c.get(2, 0) / n if n else float("nan")
        print(f"{label:<22} n={n:<8} k=2: {k2*100:5.1f}%  k>=3: {(1-k2)*100:5.1f}%  "
              f"max k={max(c) if c else 0}")
        return c

    report("snarls in A", [k for _, _, k in ba])
    report("snarls in B", [k for _, _, k in bb])
    c_pairs = report("paired tasks (A side)", [ka for ka, _ in paired])
    both2 = sum(1 for ka, kb in paired if ka == 2 and kb == 2)
    print(f"paired tasks with k=2 on both sides: {both2} / {len(paired)} = "
          f"{100*both2/len(paired) if paired else float('nan'):.1f}%")
    if a.out:
        with open(a.out, "w", newline="") as f:
            w = csv.writer(f); w.writerow(["k", "count"])
            for k in sorted(c_pairs):
                w.writerow([k, c_pairs[k]])
        print("wrote", a.out)


if __name__ == "__main__":
    main()
