#!/usr/bin/env python3
"""
graph_degree_stats.py - how many nodes of a GFA graph satisfy PanG2G's anchor
degree rule (in-degree <= 1 AND out-degree <= 1), and how long are the snarl
boundary nodes?

Reports, per graph:
  * nodes, edges, fraction of backbone (anchor-eligible) nodes,
  * length distribution of backbone nodes (how many reach the minimizer k=15),
  * length distribution of bubble boundary nodes (out-degree >= 2 sources and
    their sinks), i.e. how often boundary-sequence keys are short enough to
    collide,
  * number of distinct (source,sink) boundary keys vs. number of bubbles.

Usage:  python3 experiments/analysis/graph_degree_stats.py graph.gfa [graph2.gfa ...]
"""
import sys
from collections import Counter, defaultdict


def stats(path, kmer=15):
    seq, outd, ind = {}, Counter(), Counter()
    succ = defaultdict(list)
    n_edges = 0
    with open(path) as f:
        for line in f:
            if line.startswith("S\t"):
                p = line.rstrip("\n").split("\t")
                seq[p[1]] = p[2]
            elif line.startswith("L\t"):
                p = line.rstrip("\n").split("\t")
                outd[p[1]] += 1; ind[p[3]] += 1
                succ[p[1]].append(p[3]); n_edges += 1
    n = len(seq)
    backbone = [v for v in seq if outd[v] <= 1 and ind[v] <= 1]
    bb_len = Counter(len(seq[v]) for v in backbone)
    bb_ge_k = sum(1 for v in backbone if len(seq[v]) >= kmer)

    # bubble sources and their sinks (first common successor-of-successors)
    sources = [v for v in seq if outd[v] >= 2]
    sinks = []
    keys = Counter()
    for s in sources:
        cands = Counter()
        for a in succ[s]:
            nxt = succ.get(a, [])
            if len(nxt) == 1:
                cands[nxt[0]] += 1
            if a in cands or a in [x for x in succ[s]]:
                pass
        for a in succ[s]:
            if a in cands:
                cands[a] += 1
        if cands:
            t = cands.most_common(1)[0][0]
            sinks.append(t)
            keys[(seq[s], seq[t])] += 1
    blen = Counter(len(seq[v]) for v in sources + sinks)
    short = sum(c for L, c in blen.items() if L <= 2)
    tot_b = sum(blen.values())

    print(f"== {path}")
    print(f"  nodes={n}  edges={n_edges}  bubbles(sources with out-degree>=2)={len(sources)}")
    print(f"  anchor-eligible backbone nodes (in,out<=1): {len(backbone)} "
          f"({100*len(backbone)/n:.1f}% of nodes); with length>={kmer}: {bb_ge_k} "
          f"({100*bb_ge_k/max(1,len(backbone)):.1f}% of backbone)")
    top = sorted(bb_len.items())
    print(f"  backbone node length histogram (len:count, first 12): "
          + " ".join(f"{L}:{c}" for L, c in top[:12]) + (" ..." if len(top) > 12 else ""))
    print(f"  boundary nodes (bubble sources+sinks): {tot_b}; length<=2 bp: {short} "
          f"({100*short/max(1,tot_b):.1f}%); median len="
          f"{sorted([L for L, c in blen.items() for _ in range(c)])[tot_b//2] if tot_b else 'NA'}")
    dup = sum(c for c in keys.values() if c > 1)
    print(f"  distinct (source,sink) boundary keys: {len(keys)} for {sum(keys.values())} bubbles; "
          f"bubbles sharing a key with another bubble: {dup} "
          f"({100*dup/max(1,sum(keys.values())):.1f}%)")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        stats(p)
