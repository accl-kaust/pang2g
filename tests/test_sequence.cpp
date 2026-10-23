/**
 * tests/test_sequence.cpp
 * =======================
 * Unit tests for sequence.cpp - verifies parity with Python prototype.
 *
 * Ground truth values computed from:
 *   from g2g.sequence import edit_distance, seq_similarity
 * (run python -c "..." once and paste results here).
 */

#include <g2g/sequence.hpp>
#include <gtest/gtest.h>

using namespace g2g;

// -- edit_distance ---------------------------------------------

TEST(EditDistance, IdenticalStrings) {
    EXPECT_EQ(edit_distance("ACGT", "ACGT"), 0);
    EXPECT_EQ(edit_distance("", ""), 0);
}

TEST(EditDistance, EmptyVsNonEmpty) {
    EXPECT_EQ(edit_distance("", "ACGT"), 4);
    EXPECT_EQ(edit_distance("ACGT", ""), 4);
}

TEST(EditDistance, SingleSubstitution) {
    // ACGT → ATGT: one substitution (C→T)
    EXPECT_EQ(edit_distance("ACGT", "ATGT"), 1);
}

TEST(EditDistance, Insertion) {
    EXPECT_EQ(edit_distance("GA", "GATTACA"), 5);
}

TEST(EditDistance, DemoGraphAlleles) {
    // Values verified against Python prototype:
    //   edit_distance("GATTACA", "GATTAC") == 1  (1-bp deletion)
    //   edit_distance("CTTAGGAATCG", "GATTACA") == 8 (diverged alleles)
    EXPECT_EQ(edit_distance("GATTACA", "GATTAC"), 1);
    EXPECT_EQ(edit_distance("CTTAGGAATCG", "GATTACA"), 8);
}

// -- seq_similarity --------------------------------------------

TEST(SeqSimilarity, BothEmpty) {
    EXPECT_FLOAT_EQ(seq_similarity("", ""), 1.0f);
}

TEST(SeqSimilarity, Identical) {
    EXPECT_FLOAT_EQ(seq_similarity("ACGT", "ACGT"), 1.0f);
}

TEST(SeqSimilarity, OneSubstitution) {
    // edit=1, max_len=4 → sim = 1 - 1/4 = 0.75
    EXPECT_FLOAT_EQ(seq_similarity("ACGT", "ATGT"), 0.75f);
}

TEST(SeqSimilarity, AnchorNodes) {
    // Demo graph anchor nodes are identical → sim = 1.0
    EXPECT_FLOAT_EQ(seq_similarity("TTAGC", "TTAGC"), 1.0f);
    EXPECT_FLOAT_EQ(seq_similarity("CCCG", "CCCG"),   1.0f);
    EXPECT_FLOAT_EQ(seq_similarity("TTAAG", "TTAAG"), 1.0f);
}

TEST(SeqSimilarity, GATTACAdeletion) {
    // GATTACA vs GATTAC: edit=1, max_len=7 → sim = 1 - 1/7 ≈ 0.857
    float expected = 1.0f - 1.0f/7.0f;
    EXPECT_NEAR(seq_similarity("GATTACA", "GATTAC"), expected, 1e-5f);
}
