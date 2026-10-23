/**
 * tests/test_phase5_exp1_synthetic_gt.cpp
 * =========================================
 * Phase 5, Experiment 1 - Synthetic Ground-Truth Correctness
 *
 * Protocol (roadmap §2.7.1 / §5.1 Option A)
 * ------------------------------------------
 * 1. Build a graph pair (GA, GB) with *planted* variants at known positions:
 *      - n_shared alleles per snarl: identical sequence in A and B  → should be matched
 *      - n_a_only alleles per snarl: unique sequence in A           → should appear in a_only_alleles
 *      - n_b_only alleles per snarl: unique sequence in B           → should appear in b_only_alleles
 *
 * 2. Run the C++ ParallelAligner (Phase 4 system).
 *
 * 3. Measure recall and precision of planted variant detection:
 *      Recall_A    = |true a_only recovered| / |planted a_only|  >= 0.8
 *      Precision_A = |true a_only recovered| / |reported a_only| >= 0.8
 *      (and symmetrically for B)
 *
 * 4. Verify score agreement between ParallelAligner and sequential GraphAligner
 *    on the same inputs (rel error < 1e-3).
 *
 * Planted-graph topology
 * ----------------------
 *
 *   anc_0 ──┬── s0_sh0_a ──┬── anc_1 ──┬── s1_sh0_a ──┬── anc_2 ...
 *            ├── s0_sh1_a ──┤            ├── s1_sh1_a ──┤
 *            └── s0_ao0    ─┘            └── s1_ao0   ──┘
 *
 *   GB mirrors GA but with s{i}_sh{j}_b (same seq as s{i}_sh{j}_a)
 *   and s{i}_bo{j} (unique B-only sequences) instead of ao nodes.
 *
 * Each inner node sequence is a 12-character random DNA string seeded
 * deterministically so tests are reproducible.
 */

#include <g2g/graph.hpp>
#include <g2g/align.hpp>
#include <g2g/decompose.hpp>
#include <g2g/parallel.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <random>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

using namespace g2g;

// ── Deterministic DNA sequence generator ─────────────────────

static std::string random_dna(int length, std::mt19937& rng) {
    static const char BASES[] = "ACGT";
    std::uniform_int_distribution<int> d(0, 3);
    std::string s(length, 'A');
    for (char& c : s) c = BASES[d(rng)];
    return s;
}

// ── Planted graph description ─────────────────────────────────

/**
 * Carries ground-truth labels so tests can compute exact recall/precision.
 *
 * a_only_seqs / b_only_seqs: the inner-node sequences that were planted as
 * population-specific. These are what the aligner should report in
 * a_only_alleles / b_only_alleles respectively.
 *
 * shared_seqs: sequences present in both graphs under different node names
 * (same bytes). These should appear as matched alleles, NOT in the delta.
 */
struct PlantedGraphs {
    VariationGraph ga;
    VariationGraph gb;
    std::vector<std::string> a_only_seqs;   // planted A-only inner sequences
    std::vector<std::string> b_only_seqs;   // planted B-only inner sequences
    std::vector<std::string> shared_seqs;   // sequences that exist in both
    int n_snarls{0};
};

/**
 * Build a PlantedGraphs with:
 *   n_snarls snarls in a linear backbone chain
 *   n_shared shared alleles per snarl
 *   n_a_only A-only alleles per snarl
 *   n_b_only B-only alleles per snarl
 *   seq_len  inner-node sequence length
 *   seed     RNG seed for reproducibility
 */
static PlantedGraphs build_planted(
    int n_snarls,
    int n_shared,
    int n_a_only,
    int n_b_only,
    int seq_len = 12,
    uint32_t seed = 42)
{
    PlantedGraphs pg;
    pg.n_snarls = n_snarls;
    std::mt19937 rng(seed);

    // ── Backbone anchors - identical in both graphs ────────────
    std::vector<NodeId> anc_a, anc_b;
    for (int i = 0; i <= n_snarls; ++i) {
        std::string seq = random_dna(8, rng);
        std::string name = "anc_" + std::to_string(i);
        anc_a.push_back(pg.ga.add_node(name, seq));
        anc_b.push_back(pg.gb.add_node(name, seq));
    }

    // ── Snarl bubbles ─────────────────────────────────────────
    for (int s = 0; s < n_snarls; ++s) {
        NodeId src_a = anc_a[s],  snk_a = anc_a[s + 1];
        NodeId src_b = anc_b[s],  snk_b = anc_b[s + 1];

        // Shared alleles: same sequence, different node names
        for (int j = 0; j < n_shared; ++j) {
            std::string seq = random_dna(seq_len, rng);
            pg.shared_seqs.push_back(seq);
            std::string na = "s" + std::to_string(s) + "_sh" + std::to_string(j) + "_a";
            std::string nb = "s" + std::to_string(s) + "_sh" + std::to_string(j) + "_b";
            NodeId ia = pg.ga.add_node(na, seq);
            NodeId ib = pg.gb.add_node(nb, seq);
            pg.ga.add_edge(src_a, ia); pg.ga.add_edge(ia, snk_a);
            pg.gb.add_edge(src_b, ib); pg.gb.add_edge(ib, snk_b);
        }

        // A-only alleles
        for (int j = 0; j < n_a_only; ++j) {
            std::string seq = random_dna(seq_len, rng);
            pg.a_only_seqs.push_back(seq);
            std::string na = "s" + std::to_string(s) + "_ao" + std::to_string(j);
            NodeId ia = pg.ga.add_node(na, seq);
            pg.ga.add_edge(src_a, ia); pg.ga.add_edge(ia, snk_a);
        }

        // B-only alleles
        for (int j = 0; j < n_b_only; ++j) {
            std::string seq = random_dna(seq_len, rng);
            pg.b_only_seqs.push_back(seq);
            std::string nb = "s" + std::to_string(s) + "_bo" + std::to_string(j);
            NodeId ib = pg.gb.add_node(nb, seq);
            pg.gb.add_edge(src_b, ib); pg.gb.add_edge(ib, snk_b);
        }
    }

    // ── Representative haplotype paths (first shared allele of each snarl) ─
    std::vector<NodeId> path_a = { anc_a[0] };
    std::vector<NodeId> path_b = { anc_b[0] };
    for (int s = 0; s < n_snarls; ++s) {
        // First shared allele node name
        std::string na = "s" + std::to_string(s) + "_sh0_a";
        std::string nb = "s" + std::to_string(s) + "_sh0_b";
        path_a.push_back(pg.ga.id_map.at(na));
        path_a.push_back(anc_a[s + 1]);
        path_b.push_back(pg.gb.id_map.at(nb));
        path_b.push_back(anc_b[s + 1]);
    }
    pg.ga.add_path("hap_A0", path_a);
    pg.gb.add_path("hap_B0", path_b);

    pg.ga.finalise();
    pg.gb.finalise();
    return pg;
}

// ── Precision/Recall helpers ──────────────────────────────────

struct PRResult {
    int tp{0}, fp{0}, fn{0};
    double precision() const {
        return (tp + fp) > 0 ? static_cast<double>(tp) / (tp + fp) : 1.0;
    }
    double recall() const {
        return (tp + fn) > 0 ? static_cast<double>(tp) / (tp + fn) : 1.0;
    }
    double f1() const {
        double p = precision(), r = recall();
        return (p + r) > 0 ? 2.0 * p * r / (p + r) : 0.0;
    }
};

/**
 * Compute PR of reported delta allele sequences vs. planted ground truth.
 * We match by *sequence* (same bytes), not node name, because the aligner
 * reports DeltaAllele::sequence (the concatenated inner-node DNA string).
 */
static PRResult compute_pr(
    const std::vector<DeltaAllele>& reported,
    const std::vector<std::string>& planted_seqs)
{
    std::unordered_multiset<std::string> planted_set(planted_seqs.begin(),
                                                      planted_seqs.end());
    std::unordered_multiset<std::string> reported_set;
    for (const auto& da : reported)
        reported_set.insert(da.sequence);

    // TP = planted sequences that appear in reported
    // FN = planted that do NOT appear in reported
    // FP = reported that are NOT in planted
    PRResult r;
    auto planted_copy = planted_set;
    for (const auto& seq : reported_set) {
        auto it = planted_copy.find(seq);
        if (it != planted_copy.end()) {
            ++r.tp;
            planted_copy.erase(it);
        } else {
            ++r.fp;
        }
    }
    r.fn = static_cast<int>(planted_copy.size());
    return r;
}

// ── Test fixture ──────────────────────────────────────────────

/**
 * Parameterised by n_snarls. All tests in this suite use:
 *   2 shared alleles / snarl
 *   1 A-only allele  / snarl
 *   1 B-only allele  / snarl
 * so the expected planted sizes are:
 *   a_only_seqs.size() == n_snarls
 *   b_only_seqs.size() == n_snarls
 *   shared_seqs.size() == 2 * n_snarls
 */
class SyntheticGTTest : public ::testing::TestWithParam<int> {
protected:
    void SetUp() override {
        n_snarls_ = GetParam();
        pg_       = build_planted(n_snarls_,
                                   /*n_shared=*/2,
                                   /*n_a_only=*/1,
                                   /*n_b_only=*/1,
                                   /*seq_len=*/12,
                                   /*seed=*/42);

        // Run both aligners; results used across multiple test bodies.
        AlignerParams params;
        GraphAligner   seq_aligner(params);
        seq_result_ = seq_aligner.align(pg_.ga, pg_.gb);

        ParallelAligner par_aligner(params, /*num_threads=*/2);
        par_result_ = par_aligner.align(pg_.ga, pg_.gb);
    }

    int             n_snarls_{};
    PlantedGraphs   pg_;
    AlignmentResult seq_result_;
    AlignmentResult par_result_;
};

// ── Test bodies ───────────────────────────────────────────────

// (a) The number of snarl pairs equals n_snarls (all snarls are matched).
TEST_P(SyntheticGTTest, AllSnarlesMatched) {
    EXPECT_EQ(static_cast<int>(seq_result_.snarl_alignments.size()), n_snarls_)
        << "Expected " << n_snarls_ << " snarl alignments";
}

// (b) A-only recall >= 0.8 (planted A-only variants are found).
TEST_P(SyntheticGTTest, AOlyRecallAtLeast80Pct) {
    auto pr = compute_pr(seq_result_.a_only_alleles, pg_.a_only_seqs);
    EXPECT_GE(pr.recall(), 0.8)
        << "A-only recall too low: TP=" << pr.tp << " FN=" << pr.fn;
}

// (c) A-only precision >= 0.8 (reported A-only alleles are truly A-only).
TEST_P(SyntheticGTTest, AOnlyPrecisionAtLeast80Pct) {
    auto pr = compute_pr(seq_result_.a_only_alleles, pg_.a_only_seqs);
    EXPECT_GE(pr.precision(), 0.8)
        << "A-only precision too low: TP=" << pr.tp << " FP=" << pr.fp;
}

// (d) B-only recall >= 0.8.
TEST_P(SyntheticGTTest, BOnlyRecallAtLeast80Pct) {
    auto pr = compute_pr(seq_result_.b_only_alleles, pg_.b_only_seqs);
    EXPECT_GE(pr.recall(), 0.8)
        << "B-only recall too low: TP=" << pr.tp << " FN=" << pr.fn;
}

// (e) B-only precision >= 0.8.
TEST_P(SyntheticGTTest, BOnlyPrecisionAtLeast80Pct) {
    auto pr = compute_pr(seq_result_.b_only_alleles, pg_.b_only_seqs);
    EXPECT_GE(pr.precision(), 0.8)
        << "B-only precision too low: TP=" << pr.tp << " FP=" << pr.fp;
}

// (f) Shared alleles do NOT appear in either delta (sanity check).
TEST_P(SyntheticGTTest, SharedAllelesNotInDelta) {
    std::unordered_set<std::string> shared_set(pg_.shared_seqs.begin(),
                                                pg_.shared_seqs.end());
    for (const auto& da : seq_result_.a_only_alleles) {
        EXPECT_EQ(shared_set.count(da.sequence), 0u)
            << "Shared sequence '" << da.sequence << "' wrongly in a_only";
    }
    for (const auto& db : seq_result_.b_only_alleles) {
        EXPECT_EQ(shared_set.count(db.sequence), 0u)
            << "Shared sequence '" << db.sequence << "' wrongly in b_only";
    }
}

// (g) Total score is positive (non-trivial alignment).
TEST_P(SyntheticGTTest, TotalScorePositive) {
    EXPECT_GT(seq_result_.total_score, 0.f);
}

// (h) Score agreement between sequential and parallel aligner (< 0.1%).
TEST_P(SyntheticGTTest, ParallelMatchesSequentialScore) {
    float seq_score = seq_result_.total_score;
    float par_score = par_result_.total_score;
    if (seq_score != 0.f) {
        float rel_err = std::abs(par_score - seq_score) / std::abs(seq_score);
        EXPECT_LT(rel_err, 1e-3f)
            << "Score mismatch: seq=" << seq_score << " par=" << par_score;
    } else {
        EXPECT_NEAR(par_score, 0.f, 1e-4f);
    }
}

// (i) ParallelAligner finds the same number of snarl pairs as sequential.
TEST_P(SyntheticGTTest, ParallelSameSnarlPairCount) {
    EXPECT_EQ(par_result_.snarl_alignments.size(),
              seq_result_.snarl_alignments.size());
}

// (j) a_only and b_only collectively account for all alleles in the graph
//     (no allele goes completely unreported): |a_only| + |matched| == total A alleles.
TEST_P(SyntheticGTTest, AllAllelesAccountedFor) {
    // Total A alleles across all snarl pairs = n_snarls * (n_shared + n_a_only)
    int expected_total_a = n_snarls_ * (2 + 1);  // 2 shared + 1 a_only per snarl
    int matched_a = 0;
    for (const auto& aln : seq_result_.snarl_alignments)
        matched_a += static_cast<int>(aln.allele_matches.size());
    int a_only_count = static_cast<int>(seq_result_.a_only_alleles.size());
    EXPECT_EQ(matched_a + a_only_count, expected_total_a)
        << "matched=" << matched_a << " a_only=" << a_only_count
        << " expected_total=" << expected_total_a;
}

INSTANTIATE_TEST_SUITE_P(
    PlantedSizes,
    SyntheticGTTest,
    ::testing::Values(5, 10, 20, 50),
    [](const ::testing::TestParamInfo<int>& info) {
        return "n" + std::to_string(info.param);
    });

// ── Standalone test: planted A-only perfect recall at n=10 ────

TEST(SyntheticGTPerfect, PerfectRecallAndPrecisionN10) {
    auto pg = build_planted(10, 2, 1, 1, 12, 99);
    AlignerParams params;
    GraphAligner aligner(params);
    auto result = aligner.align(pg.ga, pg.gb);

    auto pr_a = compute_pr(result.a_only_alleles, pg.a_only_seqs);
    auto pr_b = compute_pr(result.b_only_alleles, pg.b_only_seqs);

    EXPECT_DOUBLE_EQ(pr_a.precision(), 1.0) << "Expected perfect A-only precision";
    EXPECT_DOUBLE_EQ(pr_a.recall(),    1.0) << "Expected perfect A-only recall";
    EXPECT_DOUBLE_EQ(pr_b.precision(), 1.0) << "Expected perfect B-only precision";
    EXPECT_DOUBLE_EQ(pr_b.recall(),    1.0) << "Expected perfect B-only recall";
}

// ── Edge case: single snarl ────────────────────────────────────

TEST(SyntheticGTEdge, SingleSnarlBothDeltaPopulated) {
    auto pg = build_planted(1, 2, 1, 1, 12, 7);
    AlignerParams params;
    GraphAligner aligner(params);
    auto result = aligner.align(pg.ga, pg.gb);

    EXPECT_EQ(result.snarl_alignments.size(), 1u);
    // The 1 a_only node should appear
    EXPECT_FALSE(result.a_only_alleles.empty())
        << "Expected a_only_alleles to be non-empty for single-snarl graph";
    EXPECT_FALSE(result.b_only_alleles.empty())
        << "Expected b_only_alleles to be non-empty for single-snarl graph";
}

// ── Edge case: all alleles shared (empty delta) ────────────────

TEST(SyntheticGTEdge, AllSharedEmptyDelta) {
    // n_a_only=0, n_b_only=0: everything should be matched, delta empty.
    auto pg = build_planted(5, 2, 0, 0, 12, 13);
    AlignerParams params;
    GraphAligner aligner(params);
    auto result = aligner.align(pg.ga, pg.gb);

    EXPECT_TRUE(result.a_only_alleles.empty())
        << "Expected empty a_only when all alleles shared";
    EXPECT_TRUE(result.b_only_alleles.empty())
        << "Expected empty b_only when all alleles shared";
}
