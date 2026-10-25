"""
report.py - I/O layer: result reporting
========================================
Formats and prints an AlignmentResult to stdout.
No computation here - all scores are read directly from the result object.

Functions
---------
print_result(result)   -> None   full human-readable report
"""

from __future__ import annotations
from g2g.types    import AlignmentResult
from g2g.sequence import seq_similarity, path_sequence


def print_result(result: AlignmentResult) -> None:
    ga, gb = result.graph_a, result.graph_b
    sep = "-" * 60

    print(f"\n{'='*60}")
    print(f"  Graph-to-Graph Alignment")
    print(f"  {ga.name}  ×  {gb.name}")
    print(f"{'='*60}")

    # -- Anchor nodes ------------------------------------------
    print(f"\n{sep}")
    print("  ANCHOR NODES  (shared backbone)")
    print(sep)
    if result.anchor_matches:
        for na, nb in result.anchor_matches:
            sim = seq_similarity(ga.seq(na), gb.seq(nb))
            print(f"  {na} [{ga.seq(na)!r}]  ↔  {nb} [{gb.seq(nb)!r}]"
                  f"   sim={sim:.2f}")
    else:
        print("  (none found)")

    # -- Snarl alignments --------------------------------------
    print(f"\n{sep}")
    print("  SNARL ALIGNMENTS")
    print(sep)
    for aln in result.snarl_alignments:
        sa, sb = aln.snarl_a, aln.snarl_b
        seqs_a = sa.allele_sequences(ga)
        seqs_b = sb.allele_sequences(gb)

        print(f"\n  Snarl {sa.id}  ↔  Snarl {sb.id}")
        print(f"    score={aln.score:.3f}  "
              f"(W_seq={aln.w_seq:.3f}, "
              f"W_topo={aln.w_topo:.3f}, "
              f"W_struct={aln.w_struct:.3f})")

        print("    Allele matches:")
        for i, j in aln.allele_matches:
            sim = seq_similarity(seqs_a[i], seqs_b[j])
            print(f"      A-allele[{i}] {seqs_a[i]!r}"
                  f"  ↔  B-allele[{j}] {seqs_b[j]!r}"
                  f"  sim={sim:.2f}")

        for idx in aln.a_only_alleles:
            print(f"      A-only   allele[{idx}] {seqs_a[idx]!r}"
                  f"  (no match in B)")
        for idx in aln.b_only_alleles:
            print(f"      B-only   allele[{idx}] {seqs_b[idx]!r}"
                  f"  (no match in A)")

    # -- Delta graphs ------------------------------------------
    print(f"\n{sep}")
    print("  DELTA GRAPH  ΔA→B  (A-only alleles)")
    print(sep)
    if result.a_only_alleles:
        for snarl_id, path in result.a_only_alleles:
            seq = path_sequence(ga.g, path)
            print(f"  snarl {snarl_id}  path={path}  seq={seq!r}")
    else:
        print("  (none - graphs are structurally equivalent here)")

    print(f"\n{sep}")
    print("  DELTA GRAPH  ΔB→A  (B-only alleles)")
    print(sep)
    if result.b_only_alleles:
        for snarl_id, path in result.b_only_alleles:
            seq = path_sequence(gb.g, path)
            print(f"  snarl {snarl_id}  path={path}  seq={seq!r}")
    else:
        print("  (none)")

    # -- Total score -------------------------------------------
    print(f"\n{sep}")
    print(f"  TOTAL ALIGNMENT SCORE:  {result.total_score:.4f}")
    print(sep)
