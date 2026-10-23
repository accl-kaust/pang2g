"""
types.py - Memory layer
=======================
Pure data structures with no computation.
All fields are set at construction or by the aligner; no methods
here perform any algorithmic work.

Contents
--------
Snarl           - a (source, sink) bubble with allele paths
SnarlAlignment  - alignment record for one matched snarl pair
AlignmentResult - full output: anchors, snarl alignments, delta graphs
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from g2g.graph import VariationGraph


@dataclass
class Snarl:
    """
    A snarl is a (source, sink) pair bounding a bubble subgraph S
    such that all paths from source to sink stay inside S.

    allele_paths  : list of node-id lists, each one path through the bubble
    children      : nested child snarls (snarl tree - populated by decomposer)
    """
    source: str
    sink: str
    allele_paths: list[list[str]] = field(default_factory=list)
    children: list["Snarl"]      = field(default_factory=list)

    @property
    def id(self) -> str:
        return f"({self.source},{self.sink})"

    def allele_sequences(self, vg: "VariationGraph") -> list[str]:
        """
        Return the DNA string for each allele, using inner nodes only
        (source and sink are shared anchors and excluded).
        """
        from g2g.sequence import path_sequence
        results = []
        for p in self.allele_paths:
            inner = p[1:-1]
            results.append(path_sequence(vg.g, inner) if inner else "")
        return results

    def __repr__(self) -> str:
        return (f"Snarl{self.id} "
                f"alleles={len(self.allele_paths)} "
                f"children={len(self.children)}")


@dataclass
class SnarlAlignment:
    """
    Alignment record for one matched snarl pair (sa from graph A,
    sb from graph B).

    allele_matches : list of (i, j) index pairs - allele i in sa
                     is matched to allele j in sb
    w_seq          : boundary sequence similarity contribution
    w_topo         : topological consistency contribution
    w_struct       : allele bipartite matching contribution
    score          : w_seq + w_topo + w_struct
    """
    snarl_a:       Snarl
    snarl_b:       Snarl
    score:         float
    allele_matches: list[tuple[int, int]]
    w_seq:         float
    w_topo:        float
    w_struct:      float

    @property
    def shared_allele_count(self) -> int:
        return len(self.allele_matches)

    @property
    def a_only_alleles(self) -> list[int]:
        """Indices of alleles in snarl_a that have no match in snarl_b."""
        matched_a = {i for i, _ in self.allele_matches}
        return [i for i in range(len(self.snarl_a.allele_paths))
                if i not in matched_a]

    @property
    def b_only_alleles(self) -> list[int]:
        """Indices of alleles in snarl_b that have no match in snarl_a."""
        matched_b = {j for _, j in self.allele_matches}
        return [j for j in range(len(self.snarl_b.allele_paths))
                if j not in matched_b]


@dataclass
class AlignmentResult:
    """
    Full output of a graph-to-graph alignment.

    anchor_matches    : backbone node pairs (na, nb) that are homologous
    snarl_alignments  : per-snarl alignment records
    total_score       : aggregate alignment score W(A)

    shared_nodes_a/b  : anchor node IDs in each graph  (delta layer)
    a_only_alleles    : (snarl_id, path) pairs present only in graph_a
    b_only_alleles    : (snarl_id, path) pairs present only in graph_b
    """
    graph_a:           "VariationGraph"
    graph_b:           "VariationGraph"
    snarl_alignments:  list[SnarlAlignment]
    anchor_matches:    list[tuple[str, str]]
    total_score:       float

    shared_nodes_a:  list[str]                    = field(default_factory=list)
    shared_nodes_b:  list[str]                    = field(default_factory=list)
    a_only_alleles:  list[tuple[str, list[str]]]  = field(default_factory=list)
    b_only_alleles:  list[tuple[str, list[str]]]  = field(default_factory=list)
