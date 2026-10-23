"""
demo.py - Entry point
=====================
Builds two toy population pangenomes and runs the full alignment pipeline.
Run with:  python demo.py

Demo scenario
-------------
Both graphs share:
  - Backbone nodes s0, m1, s_end  (identical sequences → anchors)
  - Bubble 1: SNP  ACGT | ATGT    (shared polymorphism)
  - Bubble 2: alleles GATTACA and GA  (shared)

Diverge at Bubble 2:
  - Population A has a third allele CTTAGGAATCG  (novel insertion, A-only)
  - Population B has a third allele GATTAC       (1-bp deletion, B-only)
"""

from __future__ import annotations
from g2g.graph    import VariationGraph
from g2g.sequence import path_sequence
from g2g.decompose import SnarlDecomposer
from g2g.align    import GraphAligner
from g2g.report   import print_result


def build_demo_graphs() -> tuple[VariationGraph, VariationGraph]:

    # -- Graph A -----------------------------------------------
    ga = VariationGraph("Population-A pangenome")

    ga.add_node("s0",    "TTAGC")    # start anchor
    ga.add_node("m1",    "CCCG")     # mid anchor
    ga.add_node("s_end", "TTAAG")    # end anchor

    ga.add_node("b1a", "ACGT")       # bubble 1: SNP allele C
    ga.add_node("b1b", "ATGT")       # bubble 1: SNP allele T

    ga.add_node("b2a", "GATTACA")    # bubble 2: shared
    ga.add_node("b2b", "GA")         # bubble 2: shared deletion
    ga.add_node("b2c", "CTTAGGAATCG")  # bubble 2: A-only novel insertion

    ga.add_path("hap_A1", ["s0", "b1a", "m1", "b2a", "s_end"])
    ga.add_path("hap_A2", ["s0", "b1b", "m1", "b2b", "s_end"])
    ga.add_path("hap_A3", ["s0", "b1a", "m1", "b2c", "s_end"])

    # -- Graph B -----------------------------------------------
    gb = VariationGraph("Population-B pangenome")

    gb.add_node("s0",    "TTAGC")
    gb.add_node("m1",    "CCCG")
    gb.add_node("s_end", "TTAAG")

    gb.add_node("b1a", "ACGT")
    gb.add_node("b1b", "ATGT")

    gb.add_node("b2a", "GATTACA")    # bubble 2: shared
    gb.add_node("b2b", "GA")         # bubble 2: shared deletion
    gb.add_node("b2d", "GATTAC")     # bubble 2: B-only 1-bp deletion

    gb.add_path("hap_B1", ["s0", "b1a", "m1", "b2a", "s_end"])
    gb.add_path("hap_B2", ["s0", "b1b", "m1", "b2b", "s_end"])
    gb.add_path("hap_B3", ["s0", "b1a", "m1", "b2d", "s_end"])

    return ga, gb


def main() -> None:
    print("\nBuilding demo variation graphs...")
    ga, gb = build_demo_graphs()
    print(ga.describe())
    print()
    print(gb.describe())

    print("\nRunning snarl decomposition...")
    snarls_a = SnarlDecomposer(ga).decompose()
    snarls_b = SnarlDecomposer(gb).decompose()
    print(f"  Graph A snarls: {snarls_a}")
    print(f"  Graph B snarls: {snarls_b}")

    print("\nRunning graph-to-graph alignment...")
    aligner = GraphAligner(alpha=1.0, beta=0.5, gamma=2.0, gap=0.5)
    result  = aligner.align(ga, gb)

    print_result(result)

    # -- Sanity checks -----------------------------------------
    print("\nSanity checks:")
    assert len(result.anchor_matches) >= 2, \
        "Should find at least 2 anchor pairs"
    assert len(result.snarl_alignments) == 2, \
        "Should align 2 snarl pairs"
    assert any(
        "CTTAGGAATCG" in path_sequence(ga.g, p)
        for _, p in result.a_only_alleles
    ), "CTTAGGAATCG insertion should be A-only"
    assert any(
        "GATTAC" in path_sequence(gb.g, p)
        for _, p in result.b_only_alleles
    ), "GATTAC deletion should be B-only"
    assert result.total_score > 0, \
        "Total score should be positive"
    print("  All assertions passed")


if __name__ == "__main__":
    main()
