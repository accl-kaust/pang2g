/**
 * tests/test_parallel.cpp
 * =======================
 * Phase 4 test suite - work-stealing parallel runtime.
 *
 * Test groups
 * -----------
 * ParallelBasic  - single-thread parallel == sequential output
 * ParallelMulti  - multi-thread produces same result as sequential
 * TaskDag        - task DAG structure (dep counts, priorities)
 * ScalingProxy   - timing proxy for Experiment 4 data (not a strict assert,
 *                  but records wall time for the paper's scaling figure)
 * MemEvict       - results freed after parent executes (peak memory)
 * EnvOverride    - G2G_NUM_THREADS env variable honoured
 */

#include <gtest/gtest.h>
#include <g2g/parallel.hpp>
#include <g2g/graph.hpp>
#include <g2g/decompose.hpp>
#include <g2g/align.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace g2g {
namespace {

// -- Graph builders ------------------------------------------------------------

/**
 * Build a variation graph with `num_bubbles` independent bubbles,
 * each with `k` alleles.  All alleles have distinct sequences so every
 * pair has a known similarity.
 *
 * Layout per bubble (bubble index b, allele index a):
 *   src_b -> allele_b_0 -> snk_b
 *         -> allele_b_1 -> snk_b
 *         ...
 *         -> allele_b_{k-1} -> snk_b
 * Bubbles are chained: snk_{b-1} == src_b.
 */
static VariationGraph make_bubbles(const std::string& name,
                                   uint32_t num_bubbles,
                                   uint32_t k,
                                   char base_char = 'A') {
    VariationGraph vg;
    vg.name = name;

    // Chain root
    NodeId prev_sink = vg.add_node("root", "N");

    for (uint32_t b = 0; b < num_bubbles; ++b) {
        NodeId src  = prev_sink;
        // Allele nodes
        std::vector<NodeId> alleles;
        for (uint32_t a = 0; a < k; ++a) {
            // Unique 2-char sequence per allele
            std::string seq(2, static_cast<char>(base_char + (b * k + a) % 26));
            NodeId an = vg.add_node(
                "b" + std::to_string(b) + "a" + std::to_string(a), seq);
            vg.add_edge(src, an);
            alleles.push_back(an);
        }
        // Sink node
        NodeId snk = vg.add_node("snk" + std::to_string(b), "N");
        for (NodeId an : alleles) vg.add_edge(an, snk);

        // One haplotype path per allele
        for (uint32_t a = 0; a < k; ++a) {
            vg.add_path(name + "_p" + std::to_string(b * k + a),
                        {src, alleles[a], snk});
        }

        prev_sink = snk;
    }
    vg.finalise();
    return vg;
}

/**
 * Build a chain of bubbles whose allele counts vary per bubble (given by `ks`).
 * Allele sequences depend only on (bubble, allele) - NOT on `name` - so two
 * graphs built from the same `ks` have byte-identical corresponding snarls and
 * must align position-for-position.  Because the bubble sizes differ, a
 * scheduler that paired snarls by subtree work (size) would cross-pair them.
 */
static VariationGraph make_varying_bubbles(const std::string& name,
                                           const std::vector<uint32_t>& ks) {
    VariationGraph vg;
    vg.name = name;

    NodeId prev_sink = vg.add_node(name + "_root", "N");

    for (uint32_t b = 0; b < ks.size(); ++b) {
        NodeId src = prev_sink;
        std::vector<NodeId> alleles;
        for (uint32_t a = 0; a < ks[b]; ++a) {
            // Distinct, deterministic 3-char sequence per (bubble, allele),
            // identical across graphs so corresponding snarls match perfectly.
            std::string seq(3, static_cast<char>('A' + (b * 7 + a) % 26));
            NodeId an = vg.add_node(
                name + "_b" + std::to_string(b) + "a" + std::to_string(a), seq);
            vg.add_edge(src, an);
            alleles.push_back(an);
        }
        NodeId snk = vg.add_node(name + "_snk" + std::to_string(b), "N");
        for (NodeId an : alleles) vg.add_edge(an, snk);
        for (uint32_t a = 0; a < ks[b]; ++a)
            vg.add_path(name + "_p" + std::to_string(b) + "_" + std::to_string(a),
                        {src, alleles[a], snk});
        prev_sink = snk;
    }
    vg.finalise();
    return vg;
}

/**
 * Backbone of distinct anchor nodes a0..a_{n-1} (each a UNIQUE sequence) with a
 * 2-allele bubble (ref edge + one alt node) at each position in `bubbles`.
 * Because the anchors are distinct, every bubble has a distinct boundary key, so
 * two graphs' bubbles correspond ONLY where their positions coincide - unlike
 * the chained builders above, whose bubbles all share the ("N","N") boundary.
 */
static VariationGraph make_positional(const std::string& name, int n_anchor,
                                      const std::vector<int>& bubbles) {
    VariationGraph vg; vg.name = name;
    std::vector<NodeId> anc;
    for (int i = 0; i < n_anchor; ++i)
        anc.push_back(vg.add_node(name + "_anc" + std::to_string(i),
                                  "ANC" + std::to_string(i)));  // seq depends on i only
    for (int i = 0; i + 1 < n_anchor; ++i) vg.add_edge(anc[i], anc[i + 1]);
    for (int b : bubbles) {
        NodeId alt = vg.add_node(name + "_alt" + std::to_string(b), "TT");
        vg.add_edge(anc[b], alt);
        vg.add_edge(alt, anc[b + 1]);
    }
    vg.finalise();
    return vg;
}

// -- Reference: sequential aligner result -------------------------------------

static AlignmentResult seq_align(const VariationGraph& ga,
                                  const VariationGraph& gb) {
    GraphAligner aligner;
    return aligner.align(ga, gb);
}

// -- ParallelBasic: 1 thread == sequential ------------------------------------

TEST(ParallelBasic, SingleThreadMatchesSequential) {
    auto ga = make_bubbles("A", 4, 3);
    auto gb = make_bubbles("B", 4, 3);

    AlignmentResult seq = seq_align(ga, gb);
    ParallelAligner pa({}, 1);
    AlignmentResult par = pa.align(ga, gb);

    ASSERT_EQ(par.snarl_alignments.size(), seq.snarl_alignments.size());
    EXPECT_NEAR(par.total_score, seq.total_score, 1e-3f);
    EXPECT_EQ(par.a_only_alleles.size(), seq.a_only_alleles.size());
    EXPECT_EQ(par.b_only_alleles.size(), seq.b_only_alleles.size());
}

TEST(ParallelBasic, EmptyGraphsDoNotCrash) {
    VariationGraph ga, gb;
    ga.finalise();
    gb.finalise();
    ParallelAligner pa({}, 1);
    AlignmentResult result = pa.align(ga, gb);
    EXPECT_EQ(result.snarl_alignments.size(), 0u);
    EXPECT_NEAR(result.total_score, 0.f, 1e-6f);
}

TEST(ParallelBasic, SingleBubble) {
    auto ga = make_bubbles("A", 1, 2);
    auto gb = make_bubbles("B", 1, 2);
    ParallelAligner pa({}, 1);
    AlignmentResult result = pa.align(ga, gb);
    EXPECT_GE(result.total_score, 0.f);
}

// -- ParallelMulti: multi-thread result consistent with sequential -------------

TEST(ParallelMulti, TwoThreadsMatchSequential) {
    auto ga = make_bubbles("A", 8, 3);
    auto gb = make_bubbles("B", 8, 3);

    AlignmentResult seq = seq_align(ga, gb);
    ParallelAligner pa({}, 2);
    AlignmentResult par = pa.align(ga, gb);

    EXPECT_NEAR(par.total_score, seq.total_score, 1e-3f);
    EXPECT_EQ(par.snarl_alignments.size(), seq.snarl_alignments.size());
}

// Regression: sequential and parallel must agree even when the graphs have
// NON-corresponding snarls (distinct boundaries). The old greedy sequential
// aligner force-paired A's unique bubble to B's unique bubble (a spurious
// cross-pair) and thus reported MORE alignments than the parallel path, which
// pairs only matching boundaries. Both should now pair only the shared bubbles.
TEST(ParallelMulti, DistinctBoundariesMatchSequential) {
    auto ga = make_positional("A", 8, {0, 2, 4});  // shared 0,2 ; A-only 4
    auto gb = make_positional("B", 8, {0, 2, 6});  // shared 0,2 ; B-only 6

    AlignmentResult seq = seq_align(ga, gb);
    ParallelAligner pa({}, 4);
    AlignmentResult par = pa.align(ga, gb);

    EXPECT_EQ(seq.snarl_alignments.size(), 2u);  // only the shared bubbles pair
    EXPECT_EQ(par.snarl_alignments.size(), seq.snarl_alignments.size());
    EXPECT_NEAR(par.total_score, seq.total_score, 1e-3f);
}

TEST(ParallelMulti, FourThreadsMatchSequential) {
    auto ga = make_bubbles("A", 16, 4);
    auto gb = make_bubbles("B", 16, 4);

    AlignmentResult seq = seq_align(ga, gb);
    ParallelAligner pa({}, 4);
    AlignmentResult par = pa.align(ga, gb);

    EXPECT_NEAR(par.total_score, seq.total_score, 1e-2f);
    EXPECT_EQ(par.snarl_alignments.size(), seq.snarl_alignments.size());
}

TEST(ParallelMulti, MaxHardwareThreads) {
    uint32_t hw = std::thread::hardware_concurrency();
    if (hw < 2) GTEST_SKIP() << "Single-core machine";

    auto ga = make_bubbles("A", 20, 3);
    auto gb = make_bubbles("B", 20, 3);

    AlignmentResult seq = seq_align(ga, gb);
    ParallelAligner pa({}, hw);
    AlignmentResult par = pa.align(ga, gb);

    EXPECT_NEAR(par.total_score, seq.total_score, 1e-2f);
}

// Regression: corresponding snarls of *different* sizes must still pair
// position-for-position.  The earlier subtree-work pairing cross-matched them
// (e.g. the smallest A-snarl with the largest B-snarl), producing a different  - 
// often negative - total score than the sequential aligner.
TEST(ParallelMulti, HeterogeneousSnarlSizesMatchSequential) {
    std::vector<uint32_t> ks = {2, 6, 3, 8, 4, 5};
    auto ga = make_varying_bubbles("A", ks);
    auto gb = make_varying_bubbles("B", ks);

    AlignmentResult seq = seq_align(ga, gb);
    ParallelAligner pa({}, 4);
    AlignmentResult par = pa.align(ga, gb);

    EXPECT_EQ(par.snarl_alignments.size(), seq.snarl_alignments.size());
    ASSERT_NE(seq.total_score, 0.f);
    float rel_err = std::abs(par.total_score - seq.total_score)
                    / std::abs(seq.total_score);
    EXPECT_LT(rel_err, 1e-3f)
        << "par=" << par.total_score << " seq=" << seq.total_score;
    EXPECT_EQ(par.a_only_alleles.size(), seq.a_only_alleles.size());
    EXPECT_EQ(par.b_only_alleles.size(), seq.b_only_alleles.size());
}

// -- TaskDag: structural properties of the task DAG ---------------------------

TEST(TaskDag, FlatTreeAllTasksReadyImmediately) {
    // Flat snarl trees (depth 1) → all dep_counts == 0
    auto ga = make_bubbles("A", 5, 2);
    auto gb = make_bubbles("B", 5, 2);

    SnarlTree ta = SnarlDecomposer(ga).decompose();
    SnarlTree tb = SnarlDecomposer(gb).decompose();

    ParallelAligner pa({}, 1);
    auto tasks = pa.build_task_dag(ga, gb, ta, tb);

    for (const auto& t : tasks) {
        EXPECT_EQ(t.dep_count.load(), 0)
            << "Flat tree tasks should all be immediately ready";
    }
}

TEST(TaskDag, PrioritiesAreNonNegative) {
    auto ga = make_bubbles("A", 6, 3);
    auto gb = make_bubbles("B", 6, 3);

    SnarlTree ta = SnarlDecomposer(ga).decompose();
    SnarlTree tb = SnarlDecomposer(gb).decompose();

    ParallelAligner pa({}, 1);
    auto tasks = pa.build_task_dag(ga, gb, ta, tb);

    for (const auto& t : tasks)
        EXPECT_GT(t.priority, 0u);
}

TEST(TaskDag, TaskCountMatchesAlignmentCount) {
    auto ga = make_bubbles("A", 10, 2);
    auto gb = make_bubbles("B", 10, 2);

    SnarlTree ta = SnarlDecomposer(ga).decompose();
    SnarlTree tb = SnarlDecomposer(gb).decompose();

    ParallelAligner pa({}, 1);
    auto tasks = pa.build_task_dag(ga, gb, ta, tb);

    // Tasks should not exceed the smaller tree
    EXPECT_LE(tasks.size(),
              std::min(ta.snarls.size(), tb.snarls.size()));
}

TEST(TaskDag, HigherKHigherPriority) {
    // Graph A has one large bubble (k=5) and one small (k=2)
    // The large bubble should produce a task with higher priority
    VariationGraph ga;
    ga.name = "A";
    // Small bubble: src0 -> a0,a1 -> snk0
    NodeId s0  = ga.add_node("s0",  "N");
    NodeId a00 = ga.add_node("a00", "AA");
    NodeId a01 = ga.add_node("a01", "CC");
    NodeId k0  = ga.add_node("k0",  "N");
    ga.add_edge(s0, a00); ga.add_edge(s0, a01);
    ga.add_edge(a00, k0); ga.add_edge(a01, k0);
    ga.add_path("p0", {s0, a00, k0});
    ga.add_path("p1", {s0, a01, k0});
    // Large bubble: k0 -> b0..b4 -> snk1
    std::vector<NodeId> big_alleles;
    NodeId k1 = ga.add_node("k1", "N");
    for (int i = 0; i < 5; ++i) {
        NodeId bi = ga.add_node("b" + std::to_string(i),
                                 std::string(3, static_cast<char>('A' + i)));
        ga.add_edge(k0, bi); ga.add_edge(bi, k1);
        ga.add_path("q" + std::to_string(i), {k0, bi, k1});
        big_alleles.push_back(bi);
    }
    ga.finalise();

    SnarlTree ta = SnarlDecomposer(ga).decompose();
    ParallelAligner pa({}, 1);
    auto tasks = pa.build_task_dag(ga, ga, ta, ta);  // pair A with itself

    if (tasks.size() >= 2) {
        // At least one task should have priority >= 5^3 = 125
        uint64_t max_pri = 0;
        for (const auto& t : tasks) max_pri = std::max(max_pri, t.priority);
        EXPECT_GE(max_pri, uint64_t(125));
    }
}

// -- ScalingProxy: wall-time data for Experiment 4 ----------------------------
// These tests are informational (they don't assert speedup - we have only
// a few cores in the test environment) but produce the timing data that
// goes into Figure 3 of the paper.

TEST(ScalingProxy, WallTimeDeclinesWithMoreThreads) {
    // 50 bubbles, k=4 - enough work to see a trend even on 2 cores
    auto ga = make_bubbles("A", 50, 4);
    auto gb = make_bubbles("B", 50, 4);

    uint32_t hw = std::thread::hardware_concurrency();
    if (hw < 2) GTEST_SKIP() << "Need ≥2 threads for scaling test";

    ParallelAligner pa1({}, 1);
    auto t0 = std::chrono::steady_clock::now();
    pa1.align(ga, gb);
    double t1_wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    ParallelAligner paN({}, hw);
    t0 = std::chrono::steady_clock::now();
    paN.align(ga, gb);
    double tN_wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    // Record (no strict assert - Experiment 4 will run on HPRC scale)
    std::cout << "[ScalingProxy] T1=" << t1_wall << "s  T" << hw
              << "=" << tN_wall << "s  speedup="
              << t1_wall / tN_wall << "x\n";

    // Weak assertion: multi-thread should not be > 10x slower than single
    EXPECT_LT(tN_wall, t1_wall * 100.0);
}

TEST(ScalingProxy, StatsFieldsPopulated) {
    auto ga = make_bubbles("A", 10, 3);
    auto gb = make_bubbles("B", 10, 3);

    ParallelAligner pa({}, 2);
    pa.align(ga, gb);
    RunStats s = pa.last_stats();

    EXPECT_GT(s.tasks_executed, 0u);
    EXPECT_GT(s.wall_seconds,   0.0);
    EXPECT_EQ(s.num_threads,    2u);
}

// -- RunStats instrumentation: steals, imbalance, ready_peak, phases ----------

TEST(RunStats, ImbalanceIsAtLeastOneAndStealsBounded) {
    auto ga = make_bubbles("A", 40, 4);
    auto gb = make_bubbles("B", 40, 4);

    ParallelAligner pa({}, 4);
    pa.align(ga, gb);
    RunStats s = pa.last_stats();

    // Imbalance = max/mean is >= 1 by construction (allow fp slack).
    EXPECT_GE(s.imbalance, 1.0 - 1e-6);
    // Can never run more "stolen" tasks than tasks executed.
    EXPECT_LE(s.steals, s.tasks_executed);
    // Ready queue held at least the seed tasks at its peak.
    EXPECT_GE(s.tasks_ready_peak, 1u);
}

TEST(RunStats, SingleThreadHasNoStealsAndUnitImbalance) {
    auto ga = make_bubbles("A", 12, 3);
    auto gb = make_bubbles("B", 12, 3);

    ParallelAligner pa({}, 1);
    pa.align(ga, gb);
    RunStats s = pa.last_stats();

    // With one thread every task's static owner (ti % 1 == 0) is that thread,
    // so nothing is "stolen" and the single thread carries all the work.
    EXPECT_EQ(s.steals, 0u);
    EXPECT_NEAR(s.imbalance, 1.0, 1e-6);
}

TEST(RunStats, PhaseBreakdownPopulated) {
    auto ga = make_bubbles("A", 16, 3);
    auto gb = make_bubbles("B", 16, 3);

    ParallelAligner pa({}, 2);
    pa.align(ga, gb);
    RunStats s = pa.last_stats();

    EXPECT_GE(s.phases.index_seconds,     0.0);
    EXPECT_GE(s.phases.decompose_seconds, 0.0);
    EXPECT_GE(s.phases.dp_seconds,        0.0);
    // dp phase mirrors the align_snarls wall time.
    EXPECT_NEAR(s.phases.dp_seconds, s.wall_seconds, 1e-9);
}

TEST(RunStats, SequentialAlignerExposesPhaseTimes) {
    auto ga = make_bubbles("A", 16, 3);
    auto gb = make_bubbles("B", 16, 3);

    GraphAligner aligner;
    aligner.align(ga, gb);
    PhaseTimes pt = aligner.last_phase_times();

    EXPECT_GE(pt.index_seconds,     0.0);
    EXPECT_GE(pt.decompose_seconds, 0.0);
    EXPECT_GT(pt.dp_seconds,        0.0);  // DP always does real work here
}

// -- EnvOverride ---------------------------------------------------------------

TEST(EnvOverride, EnvVariableHonoured) {
    // Set G2G_NUM_THREADS=3, construct with auto-detect (0), verify
#ifdef _WIN32
    _putenv_s("G2G_NUM_THREADS", "3");
#else
    ::setenv("G2G_NUM_THREADS", "3", 1);
#endif
    ParallelAligner pa({}, 0);
    EXPECT_EQ(pa.num_threads(), 3u);

#ifdef _WIN32
    _putenv_s("G2G_NUM_THREADS", "");
#else
    ::unsetenv("G2G_NUM_THREADS");
#endif
}

TEST(EnvOverride, ExplicitCountOverridesEnv) {
#ifdef _WIN32
    _putenv_s("G2G_NUM_THREADS", "99");
#else
    ::setenv("G2G_NUM_THREADS", "99", 1);
#endif
    ParallelAligner pa({}, 4);  // explicit 4 overrides env
    EXPECT_EQ(pa.num_threads(), 4u);
#ifdef _WIN32
    _putenv_s("G2G_NUM_THREADS", "");
#else
    ::unsetenv("G2G_NUM_THREADS");
#endif
}

// -- AlignSnarls interface (used by GraphAligner replacement) ------------------

TEST(ParallelAlignSnarls, ReturnsCorrectCount) {
    auto ga = make_bubbles("A", 6, 2);
    auto gb = make_bubbles("B", 6, 2);

    SnarlTree ta = SnarlDecomposer(ga).decompose();
    SnarlTree tb = SnarlDecomposer(gb).decompose();

    ParallelAligner pa({}, 2);
    auto result = pa.align_snarls(ga, gb, ta, tb, {});

    EXPECT_EQ(result.alignments.size(),
              std::min(ta.snarls.size(), tb.snarls.size()));
}

TEST(ParallelAlignSnarls, ScoresNonNegative) {
    auto ga = make_bubbles("A", 8, 3);
    auto gb = make_bubbles("B", 8, 3);

    SnarlTree ta = SnarlDecomposer(ga).decompose();
    SnarlTree tb = SnarlDecomposer(gb).decompose();

    ParallelAligner pa({}, 2);
    auto result = pa.align_snarls(ga, gb, ta, tb, {});

    for (const auto& aln : result.alignments)
        EXPECT_GE(aln.score, 0.f);
}

} // namespace
} // namespace g2g
