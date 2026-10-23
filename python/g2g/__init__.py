"""
g2g - graph-to-graph pangenome alignment
=========================================
Public API re-exports for convenient importing:

    from g2g import VariationGraph, GraphAligner, SnarlDecomposer
    from g2g import Snarl, SnarlAlignment, AlignmentResult
"""

from g2g.graph     import VariationGraph
from g2g.types     import Snarl, SnarlAlignment, AlignmentResult
from g2g.decompose import SnarlDecomposer
from g2g.align     import GraphAligner, bipartite_allele_match
from g2g.sequence  import edit_distance, seq_similarity, path_sequence
from g2g.report    import print_result

__all__ = [
    "VariationGraph",
    "Snarl",
    "SnarlAlignment",
    "AlignmentResult",
    "SnarlDecomposer",
    "GraphAligner",
    "bipartite_allele_match",
    "edit_distance",
    "seq_similarity",
    "path_sequence",
    "print_result",
]
