"""
align.py - Compute layer: graph-to-graph alignment
===================================================
Implements the alignment DP over paired snarl trees.

W(A) = W_seq  (boundary node sequence similarity)
      + W_topo (structural class consistency)
      + W_struct (bipartite allele matching within each snarl pair)

Functions
---------
bipartite_allele_match(seqs_a, seqs_b, ...)  -> (score, matches)
    Optimal partial bipartite matching between two allele sets.
    Enforces a minimum similarity threshold; unmatched alleles pay gap_cost.

Class
-----
GraphAligner(alpha, beta, gamma, gap, min_anchor_sim)
    .align(ga, gb) -> AlignmentResult

    Steps:
      1. _find_anchors       - backbone nodes with near-identical sequences
      2. SnarlDecomposer     - bubble detection on both graphs
      3. _align_snarls       - pair snarls and run per-pair DP
      4. _score_snarl_pair   - compute W_seq + W_topo + W_struct
      5. _build_deltas       - fill AlignmentResult delta fields
"""

from __future__ import annotations
from itertools import combinations, permutations
from typing import Optional

from g2g.graph     import VariationGraph
from g2g.types     import Snarl, SnarlAlignment, AlignmentResult
from g2g.sequence  import seq_similarity, path_sequence
from g2g.decompose import SnarlDecomposer


# --------------------------------------------------------------
# Bipartite allele matching
# --------------------------------------------------------------

def bipartite_allele_match(
    seqs_a: list[str],
    seqs_b: list[str],
    gap_cost: float = 0.5,
    min_match_sim: float = 0.75,
) -> tuple[float, list[tuple[int, int]]]:
    """
    Optimal weighted partial bipartite matching between two allele sets.

    Matching score for (i, j): seq_similarity(seqs_a[i], seqs_b[j])
        - forbidden if similarity < min_match_sim (prefer gap over bad match)

    Unmatched alleles in either set incur gap_cost each.

    Brute-force over all valid partial matchings of size 0..min(|A|,|B|).
    Tractable for k ≤ ~6 alleles per snarl (typical in real pangenomes).

    Returns
    -------
    (best_score, list of (i, j) index pairs)
    """
    na, nb = len(seqs_a), len(seqs_b)

    # Precompute similarity matrix
    sim = [
        [seq_similarity(seqs_a[i], seqs_b[j]) for j in range(nb)]
        for i in range(na)
    ]

    best_score = -1e9
    best_match: list[tuple[int, int]] = []

    for m in range(0, min(na, nb) + 1):
        # Choose m indices from the smaller side, permute m from the larger
        if na <= nb:
            small_range, large_range = range(na), range(nb)
            swap = False
        else:
            small_range, large_range = range(nb), range(na)
            swap = True

        for small_chosen in combinations(small_range, m):
            for large_perm in permutations(large_range, m):
                # Build (i, j) pairs - always (graph_a_idx, graph_b_idx)
                if not swap:
                    pairs = list(zip(small_chosen, large_perm))
                else:
                    pairs = [(large_perm[k], small_chosen[k]) for k in range(m)]

                # Reject pairs below similarity threshold
                if any(sim[i][j] < min_match_sim for i, j in pairs):
                    continue

                score = sum(sim[i][j] for i, j in pairs)
                score -= gap_cost * ((na - m) + (nb - m))  # unmatched on both sides

                if score > best_score:
                    best_score = score
                    best_match = list(pairs)

    if best_score == -1e9:
        # Nothing valid found - full gap penalty
        best_score = -gap_cost * (na + nb)
        best_match = []

    return best_score, best_match


# --------------------------------------------------------------
# GraphAligner
# --------------------------------------------------------------

class GraphAligner:
    """
    Aligns two VariationGraphs via anchor-seeded snarl-tree DP.

    Parameters
    ----------
    alpha          : weight for W_seq  (boundary sequence similarity)
    beta           : weight for W_topo (structural class match)
    gamma          : weight for W_struct (allele bipartite matching)
    gap            : gap penalty per unmatched allele
    min_anchor_sim : minimum sequence similarity to call two nodes anchors
    """

    def __init__(
        self,
        alpha:          float = 1.0,
        beta:           float = 0.5,
        gamma:          float = 2.0,
        gap:            float = 0.5,
        min_anchor_sim: float = 0.9,
    ):
        self.alpha          = alpha
        self.beta           = beta
        self.gamma          = gamma
        self.gap            = gap
        self.min_anchor_sim = min_anchor_sim

    # -- public entry point ------------------------------------

    def align(self, ga: VariationGraph, gb: VariationGraph) -> AlignmentResult:
        """Run full graph-to-graph alignment and return AlignmentResult."""

        # Step 1 - anchor detection
        anchors = self._find_anchors(ga, gb)

        # Step 2 - snarl decomposition
        snarls_a = SnarlDecomposer(ga).decompose()
        snarls_b = SnarlDecomposer(gb).decompose()

        # Steps 3 & 4 - pair snarls and score each pair
        snarl_alns = self._align_snarls(ga, gb, snarls_a, snarls_b, anchors)

        # Step 5 - aggregate total score
        total = sum(sa.score for sa in snarl_alns)
        for na, nb in anchors:
            total += self.alpha * seq_similarity(ga.seq(na), gb.seq(nb))

        result = AlignmentResult(
            graph_a           = ga,
            graph_b           = gb,
            snarl_alignments  = snarl_alns,
            anchor_matches    = anchors,
            total_score       = total,
        )
        self._build_deltas(result)
        return result

    # -- Step 1: anchor detection ------------------------------

    def _find_anchors(
        self, ga: VariationGraph, gb: VariationGraph
    ) -> list[tuple[str, str]]:
        """
        Match backbone nodes between graphs by near-identical sequence.
        Backbone = nodes with in-degree ≤ 1 AND out-degree ≤ 1 (chain nodes).
        Greedy best-first assignment; each node used at most once.
        """
        backbone_a = self._backbone_nodes(ga)
        backbone_b = self._backbone_nodes(gb)

        anchors:  list[tuple[str, str]] = []
        used_b:   set[str]              = set()

        for na in backbone_a:
            best_sim = self.min_anchor_sim - 1e-9
            best_nb: Optional[str] = None

            for nb in backbone_b:
                if nb in used_b:
                    continue
                s = seq_similarity(ga.seq(na), gb.seq(nb))
                if s > best_sim:
                    best_sim = s
                    best_nb  = nb

            if best_nb is not None and best_sim >= self.min_anchor_sim:
                anchors.append((na, best_nb))
                used_b.add(best_nb)

        return anchors

    def _backbone_nodes(self, vg: VariationGraph) -> list[str]:
        """Nodes with in-degree ≤ 1 and out-degree ≤ 1 - simple chain nodes."""
        g = vg.g
        return [n for n in g.nodes
                if g.in_degree(n) <= 1 and g.out_degree(n) <= 1]

    # -- Steps 3 & 4: snarl pairing ---------------------------

    def _align_snarls(
        self,
        ga:       VariationGraph,
        gb:       VariationGraph,
        snarls_a: list[Snarl],
        snarls_b: list[Snarl],
        anchors:  list[tuple[str, str]],
    ) -> list[SnarlAlignment]:
        """
        Greedy best-first pairing: for each snarl in A, find the
        highest-scoring unmatched snarl in B.
        Only snarls with a positive total score are retained.
        """
        used_b:  set[str]            = set()
        results: list[SnarlAlignment] = []

        for sa in snarls_a:
            best_score = 0.0           # threshold: only accept positive matches
            best_aln: Optional[SnarlAlignment] = None

            for sb in snarls_b:
                if sb.id in used_b:
                    continue
                aln = self._score_snarl_pair(ga, gb, sa, sb)
                if aln.score > best_score:
                    best_score = aln.score
                    best_aln   = aln

            if best_aln is not None:
                results.append(best_aln)
                used_b.add(best_aln.snarl_b.id)

        return results

    def _score_snarl_pair(
        self,
        ga: VariationGraph,
        gb: VariationGraph,
        sa: Snarl,
        sb: Snarl,
    ) -> SnarlAlignment:
        """
        Score one candidate snarl pair.

        W_seq    - average similarity of (source_a, source_b) and (sink_a, sink_b)
        W_topo   - 1.0 if both are simple bubbles (prototype assumes all are)
        W_struct - gamma * bipartite_allele_match score
        """
        # W_seq
        src_sim  = seq_similarity(ga.seq(sa.source), gb.seq(sb.source))
        sink_sim = seq_similarity(ga.seq(sa.sink),   gb.seq(sb.sink))
        w_seq    = self.alpha * (src_sim + sink_sim) / 2.0

        # W_topo  (all detected snarls are simple bubbles in this prototype)
        w_topo = self.beta * 1.0

        # W_struct
        seqs_a = sa.allele_sequences(ga)
        seqs_b = sb.allele_sequences(gb)
        struct_raw, allele_matches = bipartite_allele_match(
            seqs_a, seqs_b, gap_cost=self.gap
        )
        w_struct = self.gamma * struct_raw

        return SnarlAlignment(
            snarl_a       = sa,
            snarl_b       = sb,
            score         = w_seq + w_topo + w_struct,
            allele_matches = allele_matches,
            w_seq         = w_seq,
            w_topo        = w_topo,
            w_struct      = w_struct,
        )

    # -- Step 5: delta graph construction ---------------------

    def _build_deltas(self, result: AlignmentResult) -> None:
        """
        Populate result.shared_nodes_*, a_only_alleles, b_only_alleles
        from the completed snarl alignments.
        """
        ga, gb = result.graph_a, result.graph_b

        result.shared_nodes_a = [na for na, _  in result.anchor_matches]
        result.shared_nodes_b = [_  for _, nb  in result.anchor_matches]

        for aln in result.snarl_alignments:
            sa, sb = aln.snarl_a, aln.snarl_b

            for idx in aln.a_only_alleles:
                result.a_only_alleles.append((sa.id, sa.allele_paths[idx]))

            for idx in aln.b_only_alleles:
                result.b_only_alleles.append((sb.id, sb.allele_paths[idx]))
