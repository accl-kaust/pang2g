"""
graph.py - Memory layer: VariationGraph
========================================
Stores the variation graph in a networkx DiGraph.
All methods are either builders (mutate state) or accessors (read state).
No algorithmic computation - that lives in sequence.py / decompose.py.

Class
-----
VariationGraph
    .add_node(id, seq)          builder
    .add_edge(u, v)             builder
    .add_path(name, node_ids)   builder - also registers edges
    .seq(node_id)               accessor
    .nodes()                    accessor
    .successors(v)              accessor
    .predecessors(v)            accessor
    .describe()                 human-readable summary string
"""

from __future__ import annotations
import networkx as nx
from g2g.sequence import path_sequence


class VariationGraph:
    """
    Directed graph where every node carries a DNA sequence string.
    Edges represent adjacency in at least one haplotype path.

    Attributes
    ----------
    name   : str
    g      : nx.DiGraph  (node attr: seq)
    paths  : dict[str, list[str]]  named haplotype paths
    """

    def __init__(self, name: str = "G"):
        self.name  = name
        self.g: nx.DiGraph            = nx.DiGraph()
        self.paths: dict[str, list[str]] = {}

    # -- builders ----------------------------------------------

    def add_node(self, node_id: str, seq: str) -> "VariationGraph":
        self.g.add_node(node_id, seq=seq)
        return self

    def add_edge(self, u: str, v: str) -> "VariationGraph":
        self.g.add_edge(u, v)
        return self

    def add_path(self, name: str, node_ids: list[str]) -> "VariationGraph":
        """Register a named haplotype path and add all its edges."""
        self.paths[name] = node_ids
        for u, v in zip(node_ids, node_ids[1:]):
            self.add_edge(u, v)
        return self

    # -- accessors ---------------------------------------------

    def seq(self, node_id: str) -> str:
        return self.g.nodes[node_id]["seq"]

    def nodes(self) -> list[str]:
        return list(self.g.nodes)

    def successors(self, v: str) -> list[str]:
        return list(self.g.successors(v))

    def predecessors(self, v: str) -> list[str]:
        return list(self.g.predecessors(v))

    # -- display -----------------------------------------------

    def __repr__(self) -> str:
        return (f"VariationGraph({self.name!r}, "
                f"nodes={len(self.g)}, edges={self.g.number_of_edges()}, "
                f"paths={list(self.paths)})")

    def describe(self) -> str:
        """Human-readable dump: nodes in topological order, then paths."""
        lines = [f"=== {self.name} ==="]
        for n in nx.topological_sort(self.g):
            succs = self.successors(n)
            lines.append(
                f"  {n:6s}  seq={self.seq(n)!r:12s}  → {succs}"
            )
        for pname, pnodes in self.paths.items():
            seq = path_sequence(self.g, pnodes)
            lines.append(
                f"  path {pname!r}: {' → '.join(pnodes)}  [{seq}]"
            )
        return "\n".join(lines)
