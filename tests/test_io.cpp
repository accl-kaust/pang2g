/**
 * tests/test_io.cpp
 * =================
 * GFA reader tests — Phase 3 Task P3-8.
 *
 * Test plan (from roadmap §3.3):
 *   "Round-trip a hand-written GFA file; verify node count, edge count,
 *    and path sequences match.  Then run the demo graph as a GFA file
 *    through the full pipeline and confirm test results are unchanged."
 *
 * Suites:
 *   GfaReader   — unit tests for the parser itself.
 *   GfaRoundTrip — demo graph written as GFA, parsed back, aligner run.
 *   GfaEdgeCases — malformed input, missing sequences, W-lines.
 */

#include <g2g/io.hpp>
#include <g2g/align.hpp>
#include <g2g/sequence.hpp>
#include <gtest/gtest.h>

#include <sstream>
#include <stdexcept>
#include <string>

using namespace g2g;

// ── Helpers ───────────────────────────────────────────────────

static VariationGraph parse(const std::string& gfa_text,
                             const std::string& name = "test") {
    std::istringstream ss(gfa_text);
    return read_gfa(ss, name);
}

// ── GfaReader: basic node / edge / path parsing ───────────────

TEST(GfaReader, NodeCount) {
    auto vg = parse(
        "H\tVN:Z:1.0\n"
        "S\t1\tACGT\n"
        "S\t2\tGGGG\n"
        "S\t3\tTTTT\n"
    );
    EXPECT_EQ(vg.num_nodes(), 3u);
}

TEST(GfaReader, NodeSequences) {
    auto vg = parse(
        "S\tA\tACGT\n"
        "S\tB\tGATTACA\n"
    );
    EXPECT_EQ(vg.seq(vg.id_map.at("A")), "ACGT");
    EXPECT_EQ(vg.seq(vg.id_map.at("B")), "GATTACA");
}

TEST(GfaReader, EdgeCount) {
    auto vg = parse(
        "S\t1\tACGT\n"
        "S\t2\tGGGG\n"
        "S\t3\tTTTT\n"
        "L\t1\t+\t2\t+\t0M\n"
        "L\t2\t+\t3\t+\t0M\n"
    );
    // Edge 1→2 and 2→3
    EXPECT_EQ(vg.out_degree(vg.id_map.at("1")), 1u);
    EXPECT_EQ(vg.out_degree(vg.id_map.at("2")), 1u);
    EXPECT_EQ(vg.out_degree(vg.id_map.at("3")), 0u);
    EXPECT_EQ(vg.in_degree(vg.id_map.at("3")),  1u);
}

TEST(GfaReader, PathNodeSequence) {
    auto vg = parse(
        "S\ts\tAA\n"
        "S\tm\tCC\n"
        "S\te\tGG\n"
        "L\ts\t+\tm\t+\t0M\n"
        "L\tm\t+\te\t+\t0M\n"
        "P\thap1\ts+,m+,e+\t0M,0M\n"
    );
    ASSERT_EQ(vg.num_paths(), 1u);
    auto nodes = vg.path_nodes(0);
    ASSERT_EQ(nodes.size(), 3u);
    // Reconstruct path sequence
    std::string seq;
    for (NodeId n : nodes) seq += std::string(vg.seq(n));
    EXPECT_EQ(seq, "AACCGG");
}

TEST(GfaReader, MultiplePaths) {
    auto vg = parse(
        "S\ts0\tTTAGC\n"
        "S\tb1\tACGT\n"
        "S\tb2\tATGT\n"
        "S\te\tTTAAG\n"
        "L\ts0\t+\tb1\t+\t0M\n"
        "L\ts0\t+\tb2\t+\t0M\n"
        "L\tb1\t+\te\t+\t0M\n"
        "L\tb2\t+\te\t+\t0M\n"
        "P\thap1\ts0+,b1+,e+\t0M,0M\n"
        "P\thap2\ts0+,b2+,e+\t0M,0M\n"
    );
    EXPECT_EQ(vg.num_paths(), 2u);
    EXPECT_EQ(vg.num_nodes(), 4u);
}

TEST(GfaReader, StarSequenceStoredAsEmpty) {
    // GFA allows '*' for absent sequence
    auto vg = parse("S\tx\t*\n");
    EXPECT_EQ(vg.num_nodes(), 1u);
    EXPECT_EQ(vg.seq(0).size(), 0u);
}

TEST(GfaReader, HeaderLineIgnored) {
    auto vg = parse(
        "H\tVN:Z:1.1\txx:Z:extra\n"
        "S\tA\tACGT\n"
    );
    EXPECT_EQ(vg.num_nodes(), 1u);
}

TEST(GfaReader, CommentLineIgnored) {
    auto vg = parse(
        "# this is a comment\n"
        "S\tA\tACGT\n"
        "# another comment\n"
        "S\tB\tTTTT\n"
    );
    EXPECT_EQ(vg.num_nodes(), 2u);
}

TEST(GfaReader, OptionalTagsOnSLine) {
    // Extra tab-separated tags after the sequence field must not confuse parser
    auto vg = parse("S\tnode1\tACGT\tLN:i:4\tRC:i:10\n");
    EXPECT_EQ(vg.num_nodes(), 1u);
    EXPECT_EQ(vg.seq(0), "ACGT");
}

TEST(GfaReader, NodeNamePreservedInIdMap) {
    auto vg = parse(
        "S\tmy_node_42\tGATTACA\n"
        "S\tanother\tCCCC\n"
    );
    EXPECT_NE(vg.id_map.find("my_node_42"), vg.id_map.end());
    EXPECT_NE(vg.id_map.find("another"),    vg.id_map.end());
}

TEST(GfaReader, GraphNameSetFromArgument) {
    std::istringstream ss("S\tA\tACGT\n");
    auto vg = read_gfa(ss, "my_pangenome");
    EXPECT_EQ(vg.name, "my_pangenome");
}

// ── Reverse-complement orientation ───────────────────────────

TEST(GfaReader, NegativeOrientInPath_RcNodeCreated) {
    // A path with a '-' oriented segment should create a rc node
    auto vg = parse(
        "S\ts\tAAAA\n"
        "S\tm\tACGT\n"   // rc = ACGT (palindrome)
        "S\te\tTTTT\n"
        "L\ts\t+\tm\t-\t0M\n"
        "L\tm\t-\te\t+\t0M\n"
        "P\thap\ts+,m-,e+\t0M,0M\n"
    );
    // m_rc node should exist
    EXPECT_NE(vg.id_map.find("m_rc"), vg.id_map.end());
    NodeId m_rc = vg.id_map.at("m_rc");
    // rc of ACGT = ACGT (it's a palindrome)
    EXPECT_EQ(vg.seq(m_rc), "ACGT");
}

TEST(GfaReader, NegativeOrientInPath_RcSequenceCorrect) {
    auto vg = parse(
        "S\ta\tGATTACA\n"
        "P\tp\ta-\t*\n"
    );
    // rc of GATTACA = TGTAATC
    EXPECT_NE(vg.id_map.find("a_rc"), vg.id_map.end());
    EXPECT_EQ(vg.seq(vg.id_map.at("a_rc")), "TGTAATC");
}

// ── W-line (GFA 1.1 walk) ────────────────────────────────────

TEST(GfaReader, WalkLine_ParsedAsPath) {
    auto vg = parse(
        "S\ts1\tAAAA\n"
        "S\ts2\tCCCC\n"
        "S\ts3\tGGGG\n"
        "W\tNA12878\t1\tchr1\t0\t12\t>s1>s2>s3\n"
    );
    EXPECT_EQ(vg.num_paths(), 1u);
    // Path name format: sample#hap_idx#seq_id:start-end
    EXPECT_EQ(vg.path_names[0], "NA12878#1#chr1:0-12");
}

TEST(GfaReader, WalkLine_ReverseOrientation) {
    auto vg = parse(
        "S\ts1\tGATTACA\n"
        "S\ts2\tACGT\n"
        "W\tNA\t0\tchr1\t0\t11\t>s1<s2\n"
    );
    // s2_rc should be created for the '<' (reverse) segment
    EXPECT_NE(vg.id_map.find("s2_rc"), vg.id_map.end());
    EXPECT_EQ(vg.num_paths(), 1u);
}

// ── Error handling ────────────────────────────────────────────

TEST(GfaReader, MissingSLineFieldThrows) {
    EXPECT_THROW(parse("S\tonlyone\n"), std::runtime_error);
}

TEST(GfaReader, UnknownNodeInLLineThrows) {
    EXPECT_THROW(parse(
        "S\tA\tACGT\n"
        "L\tA\t+\tZ\t+\t0M\n"  // Z not defined
    ), std::runtime_error);
}

TEST(GfaReader, UnknownSegmentInPLineThrows) {
    EXPECT_THROW(parse(
        "S\tA\tACGT\n"
        "P\thap\tA+,MISSING+\t0M\n"
    ), std::runtime_error);
}

// ── Demo graph round-trip ─────────────────────────────────────
//
// Write the demo graphs (from test_align.cpp) as GFA strings,
// parse them back, run the aligner, and verify results match the
// hardcoded-builder versions exactly.

static const char* kDemoGfaA = R"(
H	VN:Z:1.0
S	s0	TTAGC
S	m1	CCCG
S	s_end	TTAAG
S	b1a	ACGT
S	b1b	ATGT
S	b2a	GATTACA
S	b2b	GA
S	b2c	CTTAGGAATCG
L	s0	+	b1a	+	0M
L	s0	+	b1b	+	0M
L	b1a	+	m1	+	0M
L	b1b	+	m1	+	0M
L	m1	+	b2a	+	0M
L	m1	+	b2b	+	0M
L	m1	+	b2c	+	0M
L	b2a	+	s_end	+	0M
L	b2b	+	s_end	+	0M
L	b2c	+	s_end	+	0M
P	hap_A1	s0+,b1a+,m1+,b2a+,s_end+	0M,0M,0M,0M
P	hap_A2	s0+,b1b+,m1+,b2b+,s_end+	0M,0M,0M,0M
P	hap_A3	s0+,b1a+,m1+,b2c+,s_end+	0M,0M,0M,0M
)";

static const char* kDemoGfaB = R"(
H	VN:Z:1.0
S	s0	TTAGC
S	m1	CCCG
S	s_end	TTAAG
S	b1a	ACGT
S	b1b	ATGT
S	b2a	GATTACA
S	b2b	GA
S	b2d	GATTAC
L	s0	+	b1a	+	0M
L	s0	+	b1b	+	0M
L	b1a	+	m1	+	0M
L	b1b	+	m1	+	0M
L	m1	+	b2a	+	0M
L	m1	+	b2b	+	0M
L	m1	+	b2d	+	0M
L	b2a	+	s_end	+	0M
L	b2b	+	s_end	+	0M
L	b2d	+	s_end	+	0M
P	hap_B1	s0+,b1a+,m1+,b2a+,s_end+	0M,0M,0M,0M
P	hap_B2	s0+,b1b+,m1+,b2b+,s_end+	0M,0M,0M,0M
P	hap_B3	s0+,b1a+,m1+,b2d+,s_end+	0M,0M,0M,0M
)";

class GfaRoundTrip : public ::testing::Test {
protected:
    void SetUp() override {
        std::istringstream sa(kDemoGfaA), sb(kDemoGfaB);
        ga_ = read_gfa(sa, "PopA");
        gb_ = read_gfa(sb, "PopB");
        result_ = aligner_.align(ga_, gb_);
    }
    VariationGraph ga_, gb_;
    GraphAligner   aligner_;
    AlignmentResult result_;
};

TEST_F(GfaRoundTrip, NodeCountsMatch) {
    EXPECT_EQ(ga_.num_nodes(), 8u);
    EXPECT_EQ(gb_.num_nodes(), 8u);
}

TEST_F(GfaRoundTrip, PathCountsMatch) {
    EXPECT_EQ(ga_.num_paths(), 3u);
    EXPECT_EQ(gb_.num_paths(), 3u);
}

TEST_F(GfaRoundTrip, SequencesPreserved) {
    EXPECT_EQ(ga_.seq(ga_.id_map.at("s0")),    "TTAGC");
    EXPECT_EQ(ga_.seq(ga_.id_map.at("b2c")),   "CTTAGGAATCG");
    EXPECT_EQ(gb_.seq(gb_.id_map.at("b2d")),   "GATTAC");
}

TEST_F(GfaRoundTrip, AnchorCountMatchesHardcoded) {
    // Same requirement as DemoAlignTest::FindsAtLeastTwoAnchors
    EXPECT_GE(result_.anchors.size(), 2u);
}

TEST_F(GfaRoundTrip, SnarlCountMatchesHardcoded) {
    // Same requirement as DemoAlignTest::AlignsTwoSnarlPairs
    EXPECT_EQ(result_.snarl_alignments.size(), 2u);
}

TEST_F(GfaRoundTrip, TotalScorePositive) {
    EXPECT_GT(result_.total_score, 0.f);
}

// ── write_alignment_gaf smoke test ───────────────────────────

TEST_F(GfaRoundTrip, WriteGafDoesNotThrow) {
    SnarlTree ta = SnarlDecomposer(ga_).decompose();
    SnarlTree tb = SnarlDecomposer(gb_).decompose();
    std::ostringstream out;
    EXPECT_NO_THROW(write_alignment_gaf(out, result_, ga_, gb_, ta, tb));
    // Output should contain the header line
    EXPECT_NE(out.str().find("G2G alignment GAF"), std::string::npos);
}

TEST_F(GfaRoundTrip, WriteDeltaVcfEmitsUniqueSequences) {
    std::ostringstream out;
    EXPECT_NO_THROW(write_delta_vcf(out, result_));
    const std::string s = out.str();

    // Header present
    EXPECT_NE(s.find("##fileformat=G2GdeltaVCF"), std::string::npos);

    // The ACTUAL differing sequences must appear (not a placeholder):
    //   b2c = "CTTAGGAATCG" is unique to A, b2d = "GATTAC" is unique to B.
    EXPECT_NE(s.find("CTTAGGAATCG"), std::string::npos)
        << "A-only sequence missing from delta output";
    EXPECT_NE(s.find("GATTAC"), std::string::npos)
        << "B-only sequence missing from delta output";

    // Both graphs are represented in the delta.
    EXPECT_NE(s.find("GRAPH=A"), std::string::npos);
    EXPECT_NE(s.find("GRAPH=B"), std::string::npos);
}
