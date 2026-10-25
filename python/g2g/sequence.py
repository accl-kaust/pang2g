"""
sequence.py - Compute layer: sequence operations
=================================================
All string-level computation lives here.
No data structures are defined; inputs and outputs are plain Python types.

Functions
---------
edit_distance(a, b)       -> int    Wagner-Fischer DP
seq_similarity(a, b)      -> float  normalised similarity in [0, 1]
path_sequence(graph, path) -> str   concatenate node seqs along a path
"""

from __future__ import annotations
import networkx as nx


def edit_distance(a: str, b: str) -> int:
    """
    Standard Wagner-Fischer edit distance (insert / delete / substitute).
    Space: O(n).  Time: O(m*n).
    """
    m, n = len(a), len(b)
    dp = list(range(n + 1))
    for i in range(1, m + 1):
        prev, dp[0] = dp[0], i
        for j in range(1, n + 1):
            prev, dp[j] = dp[j], min(
                prev + (0 if a[i - 1] == b[j - 1] else 1),
                dp[j] + 1,
                dp[j - 1] + 1,
            )
    return dp[n]


def seq_similarity(a: str, b: str) -> float:
    """
    Normalised similarity in [0, 1]:
        sim = 1 - edit_distance(a, b) / max(len(a), len(b))

    Returns 1.0 for identical strings, 0.0 for completely unrelated.
    Both empty strings → 1.0.
    """
    if not a and not b:
        return 1.0
    denom = max(len(a), len(b))
    return 1.0 - edit_distance(a, b) / denom


def path_sequence(graph: nx.DiGraph, path: list[str]) -> str:
    """Concatenate the 'seq' attribute of each node along path."""
    return "".join(graph.nodes[v]["seq"] for v in path)
