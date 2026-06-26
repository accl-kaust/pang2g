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

    // Graph B: bubble 1 identical (shared G); bubble 2 alt is "C" (B-only SNP)
    VariationGraph B;
    NodeId ba = B.add_node("a", "AAAAAAAA"), brc = B.add_node("rc", "C"),
           bgc = B.add_node("gc", "G"), bm = B.add_node("m", "TTTTTTTT"),
           bra = B.add_node("ra", "A"), bcc = B.add_node("cc", "C"),
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

    // Atomic variants: bubble 1 SNP C>G is shared; bubble 2 is A>T in A and
    // A>C in B (different alleles → each graph's is unique).
    ASSERT_EQ(rd->a_only_alleles.size(), 1u);
    ASSERT_EQ(rd->b_only_alleles.size(), 1u);
    EXPECT_EQ(rd->a_only_alleles[0].sequence, "T");
    EXPECT_EQ(rd->a_only_alleles[0].ref, "A");
    EXPECT_EQ(rd->b_only_alleles[0].sequence, "C");
    EXPECT_EQ(rd->a_only_alleles[0].chrom, "ref");
    EXPECT_EQ(rd->a_only_alleles[0].pos, 17);   // 8 (a) + 1 (C) + 8 (m)
    EXPECT_EQ(rd->b_only_alleles[0].pos, 17);
    // One shared atomic variant (C>G) → one matched snarl pair.
    EXPECT_EQ(rd->snarl_alignments.size(), 1u);
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

// ── PanSN / W-line reference resolution ──────────────────────

TEST(ReferenceDiff, ResolvesPanSNReferenceBySampleAndContig) {
    // Reference path named PanSN-style, as produced from a GFA W-line.
    auto build = [](const char* alt_seq) {
        VariationGraph g;
        NodeId a = g.add_node("a", "ACGTACGT"), r = g.add_node("r", "C"),
               v = g.add_node("v", alt_seq), z = g.add_node("z", "ACGTACGT");
        for (auto e : {std::pair{a, r}, {a, v}, {r, z}, {v, z}})
            g.add_edge(e.first, e.second);
        g.add_path("GRCh38#0#chr1:0-17", {a, r, z});
        g.finalise();
        return g;
    };

    VariationGraph A = build("G");
    VariationGraph B = build("G");   // same alt → nothing unique
    SnarlTree ta = SnarlDecomposer(A).decompose();
    SnarlTree tb = SnarlDecomposer(B).decompose();

    // Resolvable by sample, by contig, and by the full mangled name.
    for (const char* q : {"GRCh38", "chr1", "GRCh38#0#chr1:0-17"}) {
        auto rd = reference_diff(A, B, ta, tb, q);
        ASSERT_TRUE(rd.has_value()) << "query '" << q << "' did not resolve";
        EXPECT_EQ(rd->snarl_alignments.size(), 1u) << "query " << q;
        EXPECT_EQ(rd->a_only_alleles.size(), 0u);
    }

    // A differing alt becomes a delta tagged with the resolved contig name.
    VariationGraph C = build("T");
    SnarlTree tc = SnarlDecomposer(C).decompose();
    auto rd2 = reference_diff(A, C, ta, tc, "GRCh38");
    ASSERT_TRUE(rd2.has_value());
    ASSERT_EQ(rd2->a_only_alleles.size(), 1u);
    EXPECT_EQ(rd2->a_only_alleles[0].sequence, "G");
    EXPECT_EQ(rd2->a_only_alleles[0].chrom, "chr1");   // contig, not mangled name
    ASSERT_EQ(rd2->b_only_alleles.size(), 1u);
    EXPECT_EQ(rd2->b_only_alleles[0].sequence, "T");
}

// ── The fix: a merged bubble matches the split bubbles ───────

TEST(ReferenceDiff, MergedBubbleMatchesSplitBubbles) {
    // Graph A: ONE bubble spanning two SNPs — ref "CAT", alt "GAG"
    //   (C>G at offset 4, T>G at offset 6).
    VariationGraph A;
    NodeId a1 = A.add_node("a1", "AAAA"), cat = A.add_node("cat", "CAT"),
           gag = A.add_node("gag", "GAG"), a2 = A.add_node("a2", "GGGG");
    for (auto e : {std::pair{a1, cat}, {a1, gag}, {cat, a2}, {gag, a2}})
        A.add_edge(e.first, e.second);
    A.add_path("ref", {a1, cat, a2});
    A.finalise();

    // Graph B: the SAME two SNPs as TWO separate bubbles around a shared base.
    VariationGraph B;
    NodeId b1 = B.add_node("b1", "AAAA"), c = B.add_node("c", "C"),
           g1 = B.add_node("g1", "G"), mid = B.add_node("mid", "A"),
           t = B.add_node("t", "T"), g2 = B.add_node("g2", "G"),
           b2 = B.add_node("b2", "GGGG");
    for (auto e : {std::pair{b1, c}, {b1, g1}, {c, mid}, {g1, mid},
                   {mid, t}, {mid, g2}, {t, b2}, {g2, b2}})
        B.add_edge(e.first, e.second);
    B.add_path("ref", {b1, c, mid, t, b2});
    B.finalise();

    SnarlTree ta = SnarlDecomposer(A).decompose();
    SnarlTree tb = SnarlDecomposer(B).decompose();
    auto rd = reference_diff(A, B, ta, tb, "ref");
    ASSERT_TRUE(rd.has_value());

    // Both decompose to the SAME atomic variants (4:C>G, 6:T>G) → all shared,
    // nothing unique — even though A has one bubble and B has two.
    EXPECT_EQ(rd->a_only_alleles.size(), 0u) << "merged vs split should match";
    EXPECT_EQ(rd->b_only_alleles.size(), 0u);
    EXPECT_FLOAT_EQ(rd->total_score, 2.0f);   // two shared SNPs
}
