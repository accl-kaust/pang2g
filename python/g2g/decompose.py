"""
decompose.py - Compute layer: snarl decomposition
==================================================
Identifies bubble snarls in a VariationGraph and builds the snarl tree.
Pure computation: reads graph topology, writes Snarl objects.

Algorithm (simple bubble / dominator approximation)
----------------------------------------------------
For each node s with out-degree ≥ 2 (in topological order):
  1. BFS from each immediate successor to collect reachable sets.
  2. The sink t is the topologically earliest node in the intersection
     of all reachable sets - the first node every path from s must visit.
  3. Enumerate all simple paths from s to t as alleles.

This covers SNPs, small indels, and simple SVs.
Nested snarls: the Snarl.children field is populated when a snarl's
interior contains another snarl (detected by the same pass).

Class
-----
SnarlDecomposer(vg)
    .decompose() -> list[Snarl]
"""

from __future__ import annotations
from typing import Optional
import networkx as nx

from g2g.graph import VariationGraph
from g2g.types import Snarl


class SnarlDecomposer:
    """
    Detects simple bubble snarls in a VariationGraph.
    Requires the graph to be a DAG (no cycles).
    """

    def __init__(self, vg: VariationGraph):
        self.vg = vg

    def decompose(self) -> list[Snarl]:
        """
        Return all top-level snarls in topological order.
        Each snarl's .children are the snarls nested inside it.
        """
        g = self.vg.g

        try:
            topo = list(nx.topological_sort(g))
        except nx.NetworkXUnfeasible:
            raise ValueError("Graph contains cycles - prototype requires a DAG.")

        all_snarls: list[Snarl] = []
        visited_sources: set[str] = set()

        for s in topo:
            if s in visited_sources:
                continue
            if g.out_degree(s) < 2:
                continue

            sink = self._find_sink(s)
            if sink is None:
                continue

            alleles = self._enumerate_paths(s, sink)
            if len(alleles) < 2:
                continue

            snarl = Snarl(source=s, sink=sink, allele_paths=alleles)

            # Attach children: snarls whose source is strictly inside this one
            interior = self._interior_nodes(s, sink, alleles)
            snarl.children = [
                child for child in all_snarls
                if child.source in interior
            ]

            all_snarls.append(snarl)
            visited_sources.add(s)

        # Return only top-level snarls (not already a child of something)
        child_ids = {c.id for sn in all_snarls for c in sn.children}
        return [sn for sn in all_snarls if sn.id not in child_ids]

    # -- internal helpers --------------------------------------

    def _find_sink(self, source: str) -> Optional[str]:
        """
        Dominator approximation: find the topologically earliest node
        reachable from ALL immediate successors of source.
        """
        g = self.vg.g
        succs = list(g.successors(source))
        if len(succs) < 2:
            return None

        reach_sets: list[set[str]] = []
        for s in succs:
            r: set[str] = {s}
            queue = [s]
            while queue:
                n = queue.pop()
                for nb in g.successors(n):
                    if nb not in r:
                        r.add(nb)
                        queue.append(nb)
            reach_sets.append(r)

        candidates = reach_sets[0].intersection(*reach_sets[1:])
        if not candidates:
            return None

        try:
            topo = list(nx.topological_sort(g))
        except nx.NetworkXUnfeasible:
            return None

        for n in topo:
            if n in candidates and n != source:
                return n
        return None

    def _enumerate_paths(self, source: str, sink: str) -> list[list[str]]:
        """All simple paths from source to sink (inclusive of both endpoints)."""
        return list(nx.all_simple_paths(self.vg.g, source, sink))

    def _interior_nodes(
        self,
        source: str,
        sink: str,
        allele_paths: list[list[str]],
    ) -> set[str]:
        """All nodes that appear inside the bubble (excluding source and sink)."""
        interior: set[str] = set()
        for path in allele_paths:
            interior.update(path[1:-1])
        return interior
