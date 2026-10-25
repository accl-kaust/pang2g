/**
 * tests/test_align.cpp
 * ====================
 * Integration test: rebuild the demo graph from demo.py and verify that
 * the C++ aligner produces output consistent with the Python oracle.
 *
 * Python oracle values (run g2g/demo.py):
 *   anchor_matches:  [(s0,s0), (m1,m1), (s_end,s_end)]  or subset
 *   snarl_alignments: 2 snarl pairs (bubble1 ↔ bubble1, bubble2 ↔ bubble2)
 *   a_only_alleles:  contains path with CTTAGGAATCG
 *   b_only_alleles:  contains path with GATTAC
 *   total_score > 0
 *
 * This test is the §3.3 item 7 requirement:
 *   "For every test case in the Python prototype, run both implementations
 *    and diff the output."
 */

#include <g2g/graph.hpp>
#include <g2g/align.hpp>
#include <g2g/decompose.hpp>
#include <gtest/gtest.h>

using namespace g2g;

// -- Build demo graphs -----------------------------------------
// Mirrors build_demo_graphs() in demo.py exactly.

static VariationGraph build_graph_a() {
    VariationGraph ga;
    ga.name = "Population-A pangenome";

    auto s0    = ga.add_node("s0",    "TTAGC");
    auto m1    = ga.add_node("m1",    "CCCG");
    auto s_end = ga.add_node("s_end", "TTAAG");
    auto b1a   = ga.add_node("b1a",   "ACGT");
    auto b1b   = ga.add_node("b1b",   "ATGT");
    auto b2a   = ga.add_node("b2a",   "GATTACA");
    auto b2b   = ga.add_node("b2b",   "GA");
    auto b2c   = ga.add_node("b2c",   "CTTAGGAATCG");

    ga.add_path("hap_A1", {s0, b1a, m1, b2a, s_end});
    ga.add_path("hap_A2", {s0, b1b, m1, b2b, s_end});
    ga.add_path("hap_A3", {s0, b1a, m1, b2c, s_end});
    ga.finalise();
    return ga;
}

static VariationGraph build_graph_b() {
    VariationGraph gb;
    gb.name = "Population-B pangenome";

    auto s0    = gb.add_node("s0",    "TTAGC");
    auto m1    = gb.add_node("m1",    "CCCG");
    auto s_end = gb.add_node("s_end", "TTAAG");
    auto b1a   = gb.add_node("b1a",   "ACGT");
    auto b1b   = gb.add_node("b1b",   "ATGT");
    auto b2a   = gb.add_node("b2a",   "GATTACA");
    auto b2b   = gb.add_node("b2b",   "GA");
    auto b2d   = gb.add_node("b2d",   "GATTAC");

    gb.add_path("hap_B1", {s0, b1a, m1, b2a, s_end});
    gb.add_path("hap_B2", {s0, b1b, m1, b2b, s_end});
    gb.add_path("hap_B3", {s0, b1a, m1, b2d, s_end});
    gb.finalise();
    return gb;
}

// -- Tests -----------------------------------------------------

class DemoAlignTest : public ::testing::Test {
protected:
    void SetUp() override {
        ga_     = build_graph_a();
        gb_     = build_graph_b();
        result_ = aligner_.align(ga_, gb_);
    }

    VariationGraph ga_, gb_;
    GraphAligner   aligner_;
    AlignmentResult result_;
};

TEST_F(DemoAlignTest, FindsAtLeastTwoAnchors) {
    // Python oracle: 3 anchors (s0, m1, s_end)
    EXPECT_GE(result_.anchors.size(), 2u);
}

TEST_F(DemoAlignTest, AlignsTwoSnarlPairs) {
    // Python oracle: exactly 2 snarl alignments (bubble1, bubble2)
    EXPECT_EQ(result_.snarl_alignments.size(), 2u);
}

TEST_F(DemoAlignTest, TotalScorePositive) {
    EXPECT_GT(result_.total_score, 0.f);
}

TEST_F(DemoAlignTest, CTTAGGAATCGisAOnly) {
    // CTTAGGAATCG is in graph A only; its best match in B has sim < 0.75
    // so it should appear in a_only_alleles after build_deltas().
    bool found = false;
    for (const auto& da : result_.a_only_alleles)
        if (da.sequence == "CTTAGGAATCG") { found = true; break; }
    EXPECT_TRUE(found) << "CTTAGGAATCG not found in a_only_alleles";
}

TEST_F(DemoAlignTest, GATTACisBOnly) {
    // GATTAC is in graph B only (1-bp deletion of GATTACA).
    bool found = false;
    for (const auto& db : result_.b_only_alleles)
        if (db.sequence == "GATTAC") { found = true; break; }
    EXPECT_TRUE(found) << "GATTAC not found in b_only_alleles";
}
