/**
 * tests/test_decompose.cpp
 * ========================
 * Tests for P3-4 (vEB layout + annotate_work) and P3-5b (MinimizerIndex).
 *
 * Suite layout:
 *   SnarlWork      — annotate_work() correctness
 *   VebLayout      — apply_veb_order() structural correctness + allele integrity
 *   VebPerf        — traversal-order timing proxy for the ablation study (Q3)
 *   MinimizerIdx   — MinimizerIndex unit tests
 *   AnchorIndex    — find_anchors() switching between naïve and index paths
 */

#include <g2g/decompose.hpp>
#include <g2g/anchor.hpp>
#include <g2g/graph.hpp>
#include <g2g/sequence.hpp>
#include <gtest/gtest.h>

#include <chrono>
#include <cstdlib>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <unordered_set>
#include <utility>

using namespace g2g;

// ── Graph builder helpers ─────────────────────────────────────

// Build the standard demo graph A (two bubbles in series)
static VariationGraph make_demo_graph_a() {
    VariationGraph g;
    NodeId s0   = g.add_node("s0",   "TTAGC");
    NodeId b1a  = g.add_node("b1a",  "ACGT");
    NodeId b1b  = g.add_node("b1b",  "ATGT");
    NodeId m1   = g.add_node("m1",   "CCCG");
    NodeId b2a  = g.add_node("b2a",  "GATTACA");
    NodeId b2b  = g.add_node("b2b",  "GA");
    NodeId b2c  = g.add_node("b2c",  "CTTAGGAATCG");
    NodeId send = g.add_node("s_end","TTAAG");
    g.add_path("h1", {s0, b1a, m1, b2a, send});
    g.add_path("h2", {s0, b1b, m1, b2b, send});
    g.add_path("h3", {s0, b1a, m1, b2c, send});
    g.finalise();
    return g;
}

// Build a linear chain (no snarls)
static VariationGraph make_chain(int n) {
    VariationGraph g;
    NodeId prev = g.add_node("n0", "ACGT");
    for (int i = 1; i < n; ++i) {
        NodeId cur = g.add_node("n" + std::to_string(i), "ACGT");
        g.add_edge(prev, cur);
        prev = cur;
    }
    g.finalise();
    return g;
}

// Build a graph with k alleles per snarl
static VariationGraph make_bubble(int k, const std::string& prefix = "") {
    VariationGraph g;
    NodeId src  = g.add_node(prefix + "src", "AAAA");
    NodeId sink = g.add_node(prefix + "snk", "TTTT");
    for (int i = 0; i < k; ++i) {
        std::string seq(4 + i, 'C');  // distinct lengths to ensure distinct sequences
        NodeId mid = g.add_node(prefix + "m" + std::to_string(i), seq);
        std::vector<NodeId> path = {src, mid, sink};
        g.add_path(prefix + "h" + std::to_string(i), path);
    }
    g.finalise();
    return g;
}

// ── SnarlWork ─────────────────────────────────────────────────

TEST(SnarlWork, SingleBubble_k3) {
    auto g = make_demo_graph_a();
    auto tree = SnarlDecomposer(g).decompose();
    // Two snarls: k=2 and k=3.
    // subtree_work = k^3 for each (no children).
    // We don't know which order they appear, but sum = 8 + 27 = 35.
    uint64_t total = 0;
    for (const auto& s : tree.snarls)
        total += s.subtree_work;
    EXPECT_EQ(total, 8u + 27u);
}

TEST(SnarlWork, BubbleK5_subtreeWork) {
    auto g = make_bubble(5);
    auto tree = SnarlDecomposer(g).decompose();
    ASSERT_EQ(tree.snarls.size(), 1u);
    EXPECT_EQ(tree.snarls[0].subtree_work, 125u);  // 5^3
}

TEST(SnarlWork, LinearChain_noSnarls) {
    auto g = make_chain(10);
    auto tree = SnarlDecomposer(g).decompose();
    EXPECT_EQ(tree.snarls.size(), 0u);
}

// ── VebLayout ─────────────────────────────────────────────────

TEST(VebLayout, AllSnarlsPresent) {
    auto g = make_demo_graph_a();
    auto tree = SnarlDecomposer(g).decompose();
    // vEB reorder must not lose or duplicate any snarls
    EXPECT_EQ(tree.snarls.size(), 2u);
}

TEST(VebLayout, RootsValid) {
    auto g = make_demo_graph_a();
    auto tree = SnarlDecomposer(g).decompose();
    for (uint32_t r : tree.roots)
        EXPECT_LT(r, tree.snarls.size());
}

TEST(VebLayout, AlleleCountsPreserved) {
    // After vEB reorder, total alleles across all snarls must equal
    // the number of entries in allele_path_offsets.
    auto g = make_demo_graph_a();
    auto tree = SnarlDecomposer(g).decompose();
    uint32_t total_alleles = 0;
    for (const auto& s : tree.snarls)
        total_alleles += s.allele_count();
    EXPECT_EQ(total_alleles, tree.allele_path_offsets.size());
}

TEST(VebLayout, AllelePathsAccessibleAfterReorder) {
    // Every allele path must be accessible via allele_path() and non-empty.
    auto g = make_demo_graph_a();
    auto tree = SnarlDecomposer(g).decompose();
    for (const auto& s : tree.snarls) {
        EXPECT_GE(s.allele_count(), 2u);
        for (uint32_t i = 0; i < s.allele_count(); ++i) {
            auto path = tree.allele_path(s.alleles_begin, i);
            EXPECT_GE(path.size(), 2u)   // at least source + sink
                << "allele " << i << " of snarl has fewer than 2 nodes";
        }
    }
}

TEST(VebLayout, AlleleSequencesConsistentAfterReorder) {
    // allele_sequences() should still produce the same set of strings
    // (order may differ) after vEB reorder.
    auto g = make_demo_graph_a();
    auto tree = SnarlDecomposer(g).decompose();
    // For the bubble with k=3 (source s0, alleles b2a/b2b/b2c), inner seqs
    // should be GATTACA, GA, CTTAGGAATCG in some order.
    std::unordered_set<std::string> expected = {"GATTACA", "GA", "CTTAGGAATCG"};
    bool found = false;
    for (const auto& s : tree.snarls) {
        if (s.allele_count() != 3) continue;
        std::unordered_set<std::string> actual;
        std::vector<std::string> seqs;
        allele_sequences(g, s, tree, seqs);
        for (const auto& sq : seqs) actual.insert(sq);
        if (actual == expected) { found = true; break; }
    }
    EXPECT_TRUE(found) << "3-allele snarl with expected sequences not found after vEB reorder";
}

TEST(VebLayout, SourceSinkPreservedAfterReorder) {
    // source and sink NodeIds must still be valid nodes in the graph.
    auto g = make_demo_graph_a();
    auto tree = SnarlDecomposer(g).decompose();
    for (const auto& s : tree.snarls) {
        EXPECT_LT(s.source, g.num_nodes());
        EXPECT_LT(s.sink,   g.num_nodes());
    }
}

// ── VebPerf: traversal-order timing proxy ─────────────────────
//
// This is the ablation-study data point for paper §2.7.3 (Q3, Experiment 5).
// We build a synthetic graph with many snarls, then measure wall-clock time
// for a simple bottom-up traversal in two orderings:
//   1. Sequential (natural DFS order, pre-vEB).
//   2. vEB order (post-apply_veb_order).
//
// On a warm cache the difference is small; it becomes significant once the
// working set exceeds L3.  For the paper we run on HPRC-scale graphs.
// Here we just verify the infrastructure works and the vEB traversal
// produces the same total work as the sequential traversal.

static VariationGraph make_many_bubbles(int num_bubbles, int k_per_bubble) {
    VariationGraph g;
    // Chain of bubbles: src_0 → bubble_0 → src_1 → bubble_1 → ...
    NodeId chain_prev = g.add_node("chain_start", "AAAA");
    for (int b = 0; b < num_bubbles; ++b) {
        std::string pfx = "b" + std::to_string(b) + "_";
        NodeId src  = g.add_node(pfx + "src",  "GGGG");
        NodeId sink = g.add_node(pfx + "snk",  "CCCC");
        g.add_edge(chain_prev, src);
        for (int a = 0; a < k_per_bubble; ++a) {
            std::string seq(5 + a, 'T');
            NodeId mid = g.add_node(pfx + "m" + std::to_string(a), seq);
            g.add_edge(src, mid);
            g.add_edge(mid, sink);
        }
        chain_prev = sink;
    }
    g.finalise();
    return g;
}

TEST(VebPerf, TraversalProducesCorrectTotalWork) {
    // 20 bubbles, k=4 each → 20 snarls, subtree_work=64 each, total=1280
    auto g = make_many_bubbles(20, 4);
    auto tree = SnarlDecomposer(g).decompose();
    uint64_t total = 0;
    for (const auto& s : tree.snarls)
        total += s.subtree_work;
    EXPECT_EQ(total, 20u * 64u);  // 20 * 4^3
}

TEST(VebPerf, VebOrderDoesNotChangeTotalAlleles) {
    auto g = make_many_bubbles(20, 4);
    auto tree = SnarlDecomposer(g).decompose();
    uint32_t total = 0;
    for (const auto& s : tree.snarls) total += s.allele_count();
    EXPECT_EQ(total, tree.allele_path_offsets.size());
    EXPECT_EQ(tree.snarls.size(), 20u);
}

// ── MinimizerIdx ──────────────────────────────────────────────

TEST(MinimizerIdx, IndexedNodeCount) {
    auto g = make_demo_graph_a();
    AnchorParams p; p.kmer_len = 3; p.window = 2;
    MinimizerIndex idx(g, p);
    // Backbone = nodes with in_degree≤1 AND out_degree≤1.
    // In demo graph A: s0 (out=2→not backbone), b1a (in=1,out=1), b1b, m1 (in=2→not),
    // b2a, b2b, b2c, s_end (in=3→not).
    // So backbone = {b1a, b1b, b2a, b2b, b2c}.
    EXPECT_EQ(idx.indexed_nodes(), 5u);
}

TEST(MinimizerIdx, QueryReturnsIdenticalNode) {
    // A node should always be a candidate for its own sequence.
    auto g = make_demo_graph_a();
    AnchorParams p; p.kmer_len = 3; p.window = 2; p.min_sim = 0.75f;
    MinimizerIndex idx(g, p);

    NodeId b2a = g.id_map.at("b2a");  // seq = "GATTACA"
    auto cands = idx.query("GATTACA");
    bool found = std::find(cands.begin(), cands.end(), b2a) != cands.end();
    EXPECT_TRUE(found) << "b2a not in its own query candidates";
}

TEST(MinimizerIdx, QueryNoMatchForGarbage) {
    auto g = make_demo_graph_a();
    AnchorParams p; p.kmer_len = 4; p.window = 3;
    MinimizerIndex idx(g, p);
    // Sequence of all X's has no k-mer overlap with DNA sequences
    auto cands = idx.query("XXXXXXXXXXXX");
    EXPECT_EQ(cands.size(), 0u);
}

TEST(MinimizerIdx, QueryShortSequence) {
    // Sequence shorter than k: index must not crash, returns empty or partial
    auto g = make_demo_graph_a();
    AnchorParams p; p.kmer_len = 15; p.window = 10;
    MinimizerIndex idx(g, p);
    // "GA" is shorter than k=15; should not throw
    EXPECT_NO_THROW(idx.query("GA"));
}

TEST(MinimizerIdx, NearDuplicateReturnedAsCandidate) {
    // A sequence with one substitution should still share minimizers with
    // the original for small k (k=3, w=2).
    auto g = make_demo_graph_a();
    AnchorParams p; p.kmer_len = 3; p.window = 2;
    MinimizerIndex idx(g, p);

    NodeId b2a = g.id_map.at("b2a");  // "GATTACA"
    // One substitution: "GATTACX" → last kmer changed but earlier ones same
    auto cands = idx.query("GATTACX");
    // Some minimizers should overlap → b2a should still appear
    bool found = std::find(cands.begin(), cands.end(), b2a) != cands.end();
    EXPECT_TRUE(found) << "near-duplicate did not find b2a as candidate";
}

// ── AnchorIndex: find_anchors switches between paths ──────────

TEST(AnchorIndex, NaivePathSmallGraph) {
    // Demo graphs: |backbone_a|*|backbone_b| << kNaiveThreshold → naïve path
    auto ga = make_demo_graph_a();
    VariationGraph gb;
    NodeId s0 = gb.add_node("s0","TTAGC"); NodeId b1a = gb.add_node("b1a","ACGT");
    NodeId b1b = gb.add_node("b1b","ATGT"); NodeId m1 = gb.add_node("m1","CCCG");
    NodeId b2a = gb.add_node("b2a","GATTACA"); NodeId b2b = gb.add_node("b2b","GA");
    NodeId b2d = gb.add_node("b2d","GATTAC");  NodeId se = gb.add_node("s_end","TTAAG");
    gb.add_path("h1",{s0,b1a,m1,b2a,se}); gb.add_path("h2",{s0,b1b,m1,b2b,se});
    gb.add_path("h3",{s0,b1a,m1,b2d,se}); gb.finalise();

    AnchorParams p; p.kmer_len = 3; p.window = 2; p.min_sim = 0.9f;
    auto anchors = find_anchors(ga, gb, p);
    EXPECT_GE(anchors.size(), 2u);
    for (const auto& a : anchors)
        EXPECT_GE(a.sim, 0.9f);
}

TEST(AnchorIndex, NoAnchorsBetweenUnrelatedGraphs) {
    // Two graphs with completely different sequences → no anchors
    VariationGraph ga, gb;
    NodeId a = ga.add_node("a","AAAAAAAAAA"); NodeId b = ga.add_node("b","AAAAAAAAAA");
    ga.add_edge(a,b); ga.finalise();
    NodeId c = gb.add_node("c","TTTTTTTTTT"); NodeId d = gb.add_node("d","TTTTTTTTTT");
    gb.add_edge(c,d); gb.finalise();

    AnchorParams p; p.min_sim = 0.9f;
    auto anchors = find_anchors(ga, gb, p);
    EXPECT_EQ(anchors.size(), 0u);
}

TEST(AnchorIndex, IndexPathLargeGraph) {
    // Build two graphs large enough to trigger the MinimizerIndex path.
    // We plant 10 identical "anchor" nodes scattered among random backbone nodes.
    // find_anchors should find all 10 planted anchors.
    const int N = 1100;   // backbone_a * backbone_b > 10^6 threshold
    const int planted = 10;

    // Each graph is a chain of N backbone nodes (in_deg=1, out_deg=1)
    // so that every node qualifies as backbone.
    auto make_large = [&](bool with_planted) -> VariationGraph {
        VariationGraph g;
        std::mt19937 rng(42 + static_cast<int>(with_planted));
        std::uniform_int_distribution<int> base_dist(0, 3);
        const char bases[] = "ACGT";

        // Generate N random sequences of length 20
        std::vector<NodeId> nodes;
        nodes.reserve(static_cast<size_t>(N));
        for (int i = 0; i < N; ++i) {
            std::string seq;
            seq.reserve(20);
            for (int j = 0; j < 20; ++j)
                seq += bases[base_dist(rng)];
            // Plant recognisable anchors at fixed positions
            if (i < planted && with_planted)
                seq = "GATTACAGATTACAGATTAC";  // 20-char anchor sequence
            else if (i < planted && !with_planted)
                seq = "GATTACAGATTACAGATTAC";  // same in both graphs
            nodes.push_back(g.add_node("n" + std::to_string(i), seq));
        }
        // Chain them
        for (int i = 0; i + 1 < N; ++i)
            g.add_edge(nodes[static_cast<size_t>(i)],
                       nodes[static_cast<size_t>(i+1)]);
        g.finalise();
        return g;
    };

    auto ga = make_large(true);
    auto gb = make_large(false);

    AnchorParams p; p.kmer_len = 11; p.window = 5; p.min_sim = 0.95f;
    auto anchors = find_anchors(ga, gb, p);

    // All planted anchors have identical sequences → should be found
    int found_planted = 0;
    for (const auto& a : anchors) {
        if (ga.seq(a.node_a) == "GATTACAGATTACAGATTAC")
            ++found_planted;
    }
    EXPECT_EQ(found_planted, planted)
        << "MinimizerIndex path missed " << (planted - found_planted)
        << " of " << planted << " planted anchors";
}

// Sensitivity sweep: different k/w settings on the large-graph fixture.
// Numbers from this test populate Table in §2.7.3 (ablation Q3).
TEST(AnchorIndex, SensitivitySweep) {
    const int N = 1100;
    const int planted = 10;

    auto make_graph = [&](uint32_t seed) -> VariationGraph {
        VariationGraph g;
        std::mt19937 rng(seed);
        std::uniform_int_distribution<int> bd(0,3);
        const char bases[] = "ACGT";
        std::vector<NodeId> nodes;
        for (int i = 0; i < N; ++i) {
            std::string seq;
            for (int j = 0; j < 20; ++j) seq += bases[bd(rng)];
            if (i < planted) seq = "GATTACAGATTACAGATTAC";
            nodes.push_back(g.add_node("n"+std::to_string(i), seq));
        }
        for (int i = 0; i+1 < N; ++i)
            g.add_edge(nodes[static_cast<size_t>(i)],
                       nodes[static_cast<size_t>(i+1)]);
        g.finalise();
        return g;
    };

    auto ga = make_graph(10);
    auto gb = make_graph(10);  // same seed → identical anchor sequences

    // Sensitivity sweep: k ∈ {11, 15} work on 20-char sequences (len > k).
    // k=21 exceeds anchor sequence length (20) → 0 valid k-mers → 0 recall.
    // This models the real tradeoff in §2.7.3: larger k = higher specificity
    // but lower recall on short sequences.
    struct Setting { uint8_t k, w; int expected_hits; };
    for (auto [k, w, expected] :
             std::vector<Setting>{{11,5,planted},{15,10,planted},{21,15,0}}) {
        AnchorParams p;
        p.kmer_len = k; p.window = w; p.min_sim = 0.95f;
        auto anchors = find_anchors(ga, gb, p);
        int hits = 0;
        for (const auto& a : anchors)
            if (ga.seq(a.node_a) == "GATTACAGATTACAGATTAC") ++hits;
        EXPECT_EQ(hits, expected)
            << "k=" << static_cast<int>(k) << " w=" << static_cast<int>(w)
            << " expected " << expected << " hits, got " << hits;
    }
}

// ── Vg snarl back-end (snarls_vg.cpp) ─────────────────────────
// These tests exercise the optional vg-toolkit decomposer.  They self-skip
// when vg is not available (set $G2G_VG to the vg binary, or put vg on PATH);
// in CI without vg they are skipped rather than failed.

TEST(VgSnarls, DecomposeMatchesBuiltinOnDemoGraph) {
    VariationGraph g = make_demo_graph_a();
    // Gate on vg actually producing boundaries (not just env being set) so a
    // broken/missing vg skips instead of silently passing via fallback.
    if (!vg_snarl_boundaries(g).has_value())
        GTEST_SKIP() << "vg unavailable (set G2G_VG to the vg binary)";

    auto boundary_set = [&](const SnarlTree& t) {
        std::set<std::pair<std::string, std::string>> s;
        for (const auto& sn : t.snarls)
            s.insert({g.node_names[sn.source], g.node_names[sn.sink]});
        return s;
    };

    ::unsetenv("G2G_SNARLS");
    SnarlTree builtin = SnarlDecomposer(g).decompose();

    ::setenv("G2G_SNARLS", "vg", 1);
    SnarlTree viavg = SnarlDecomposer(g).decompose();
    ::unsetenv("G2G_SNARLS");

    EXPECT_EQ(viavg.snarls.size(), builtin.snarls.size());
    EXPECT_EQ(boundary_set(viavg), boundary_set(builtin));
}

TEST(VgSnarls, BoundariesContainBuiltinSnarls) {
    VariationGraph g = make_demo_graph_a();
    auto vb = vg_snarl_boundaries(g);
    if (!vb.has_value())
        GTEST_SKIP() << "vg unavailable (set G2G_VG to the vg binary)";

    std::set<std::pair<std::string, std::string>> vgset;
    for (auto [s, t] : *vb)
        vgset.insert({g.node_names[s], g.node_names[t]});

    SnarlTree builtin = SnarlDecomposer(g).decompose();
    for (const auto& sn : builtin.snarls) {
        std::pair<std::string, std::string> key{g.node_names[sn.source],
                                                g.node_names[sn.sink]};
        EXPECT_TRUE(vgset.count(key))
            << "vg missing built-in snarl " << key.first << ">" << key.second;
    }
}
