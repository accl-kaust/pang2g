/**
 * tests/test_refdiff.cpp
 * ======================
 * Reference-coordinate graph-to-graph diff (refdiff.cpp).
 *
 * These tests are pure G2G (no vg dependency): they build small graphs that
 * share an embedded reference path and verify that variants are matched by
 * reference position — including when the two graphs chop the shared reference
 * into different nodes, which is exactly what breaks boundary-sequence
 * correspondence on real pangenome graphs.
 */

#include <g2g/graph.hpp>
#include <g2g/decompose.hpp>
#include <g2g/refdiff.hpp>
#include <gtest/gtest.h>

#include <string>

using namespace g2g;

// ── Shared variant + one A-only + one B-only ─────────────────

TEST(ReferenceDiff, SharedAndUniqueVariants) {
    // Graph A: a -[C|G]- m -[A|T]- z   (reference path goes C, A)
    VariationGraph A;
    NodeId a = A.add_node("a", "AAAAAAAA"), rc = A.add_node("rc", "C"),
           gc = A.add_node("gc", "G"), m = A.add_node("m", "TTTTTTTT"),
           ra = A.add_node("ra", "A"), tx = A.add_node("tx", "T"),
           z = A.add_node("z", "GGGGGGGG");
    for (auto e : {std::pair{a, rc}, {a, gc}, {rc, m}, {gc, m},
                   {m, ra}, {m, tx}, {ra, z}, {tx, z}})
        A.add_edge(e.first, e.second);
    A.add_path("ref", {a, rc, m, ra, z});
    A.finalise();

    // Graph B: bubble 1 identical (shared G); bubble 2 alt is "CC" (B-only)
    VariationGraph B;
    NodeId ba = B.add_node("a", "AAAAAAAA"), brc = B.add_node("rc", "C"),
           bgc = B.add_node("gc", "G"), bm = B.add_node("m", "TTTTTTTT"),
           bra = B.add_node("ra", "A"), bcc = B.add_node("cc", "CC"),
           bz = B.add_node("z", "GGGGGGGG");
    for (auto e : {std::pair{ba, brc}, {ba, bgc}, {brc, bm}, {bgc, bm},
                   {bm, bra}, {bm, bcc}, {bra, bz}, {bcc, bz}})
        B.add_edge(e.first, e.second);
    B.add_path("ref", {ba, brc, bm, bra, bz});
    B.finalise();

    SnarlTree ta = SnarlDecomposer(A).decompose();
    SnarlTree tb = SnarlDecomposer(B).decompose();
    auto rd = reference_diff(A, B, ta, tb, "ref");
    ASSERT_TRUE(rd.has_value());

    ASSERT_EQ(rd->a_only_alleles.size(), 1u);
    ASSERT_EQ(rd->b_only_alleles.size(), 1u);
    EXPECT_EQ(rd->a_only_alleles[0].sequence, "T");
    EXPECT_EQ(rd->b_only_alleles[0].sequence, "CC");
    EXPECT_EQ(rd->a_only_alleles[0].chrom, "ref");
    EXPECT_EQ(rd->a_only_alleles[0].pos, 17);   // 8 (a) + 1 (C) + 8 (m)
    EXPECT_EQ(rd->b_only_alleles[0].pos, 17);
    // Both sites present in both graphs → 2 alignments; one shared alt (G).
    EXPECT_EQ(rd->snarl_alignments.size(), 2u);
    EXPECT_FLOAT_EQ(rd->total_score, 1.0f);
}

// ── Robustness: shared variant matches despite different segmentation ─

TEST(ReferenceDiff, RobustToDifferentSegmentation) {
    // A: one 12bp anchor before the bubble.
    VariationGraph A;
    NodeId a = A.add_node("a", "AAAAAAAAAAAA"), rc = A.add_node("rc", "C"),
           gc = A.add_node("gc", "G"), z = A.add_node("z", "GGGGGGGGGGGG");
    for (auto e : {std::pair{a, rc}, {a, gc}, {rc, z}, {gc, z}})
        A.add_edge(e.first, e.second);
    A.add_path("ref", {a, rc, z});
    A.finalise();

    // B: the SAME leading 12bp split into 5 + 7 (different node boundaries),
    // same variant (G) at the same reference position.
    VariationGraph B;
    NodeId a1 = B.add_node("a1", "AAAAA"), a2 = B.add_node("a2", "AAAAAAA"),
           rc2 = B.add_node("rc", "C"), gc2 = B.add_node("gc", "G"),
           z2 = B.add_node("z", "GGGGGGGGGGGG");
    for (auto e : {std::pair{a1, a2}, {a2, rc2}, {a2, gc2}, {rc2, z2}, {gc2, z2}})
        B.add_edge(e.first, e.second);
    B.add_path("ref", {a1, a2, rc2, z2});
    B.finalise();

    SnarlTree ta = SnarlDecomposer(A).decompose();
    SnarlTree tb = SnarlDecomposer(B).decompose();
    auto rd = reference_diff(A, B, ta, tb, "ref");
    ASSERT_TRUE(rd.has_value());

    // The bubble opens at reference offset 12 in BOTH (A: 12; B: 5+7) → matched,
    // shared alt G, nothing unique — even though node boundaries differ.
    EXPECT_EQ(rd->a_only_alleles.size(), 0u);
    EXPECT_EQ(rd->b_only_alleles.size(), 0u);
    EXPECT_EQ(rd->snarl_alignments.size(), 1u);
    EXPECT_FLOAT_EQ(rd->total_score, 1.0f);
}

// ── Missing reference path → nullopt ─────────────────────────

TEST(ReferenceDiff, MissingReferencePathReturnsNullopt) {
    VariationGraph A;
    NodeId a = A.add_node("a", "ACGT"), b = A.add_node("b", "ACGT");
    A.add_edge(a, b);
    A.add_path("ref", {a, b});
    A.finalise();

    SnarlTree ta = SnarlDecomposer(A).decompose();
    EXPECT_FALSE(reference_diff(A, A, ta, ta, "no_such_path").has_value());
}
