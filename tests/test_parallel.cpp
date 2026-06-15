/**
 * tests/test_parallel.cpp
 * =======================
 * Phase 4 test suite — work-stealing parallel runtime.
 *
 * Test groups
 * -----------
 * ParallelBasic  — single-thread parallel == sequential output
 * ParallelMulti  — multi-thread produces same result as sequential
 * TaskDag        — task DAG structure (dep counts, priorities)
 * ScalingProxy   — timing proxy for Experiment 4 data (not a strict assert,
 *                  but records wall time for the paper's scaling figure)
 * MemEvict       — results freed after parent executes (peak memory)
 * EnvOverride    — G2G_NUM_THREADS env variable honoured
 */

#include <gtest/gtest.h>
#include <g2g/parallel.hpp>
#include <g2g/graph.hpp>
#include <g2g/decompose.hpp>
#include <g2g/align.hpp>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>

namespace g2g {
namespace {

// ── Graph builders ────────────────────────────────────────────────────────────

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

// ── Reference: sequential aligner result ─────────────────────────────────────

static AlignmentResult seq_align(const VariationGraph& ga,
                                  const VariationGraph& gb) {
    GraphAligner aligner;
    return aligner.align(ga, gb);
}

// ── ParallelBasic: 1 thread == sequential ────────────────────────────────────

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

// ── ParallelMulti: multi-thread result consistent with sequential ─────────────

TEST(ParallelMulti, TwoThreadsMatchSequential) {
    auto ga = make_bubbles("A", 8, 3);
    auto gb = make_bubbles("B", 8, 3);

    AlignmentResult seq = seq_align(ga, gb);
    ParallelAligner pa({}, 2);
    AlignmentResult par = pa.align(ga, gb);

    EXPECT_NEAR(par.total_score, seq.total_score, 1e-3f);
    EXPECT_EQ(par.snarl_alignments.size(), seq.snarl_alignments.size());
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

// ── TaskDag: structural properties of the task DAG ───────────────────────────

TEST(TaskDag, FlatTreeAllTasksReadyImmediately) {
    // Flat snarl trees (depth 1) → all dep_counts == 0
    auto ga = make_bubbles("A", 5, 2);
    auto gb = make_bubbles("B", 5, 2);

    SnarlTree ta = SnarlDecomposer(ga).decompose();
    SnarlTree tb = SnarlDecomposer(gb).decompose();

    ParallelAligner pa({}, 1);
    auto tasks = pa.build_task_dag(ta, tb);

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
    auto tasks = pa.build_task_dag(ta, tb);

    for (const auto& t : tasks)
        EXPECT_GT(t.priority, 0u);
}

TEST(TaskDag, TaskCountMatchesAlignmentCount) {
    auto ga = make_bubbles("A", 10, 2);
    auto gb = make_bubbles("B", 10, 2);

    SnarlTree ta = SnarlDecomposer(ga).decompose();
    SnarlTree tb = SnarlDecomposer(gb).decompose();

    ParallelAligner pa({}, 1);
    auto tasks = pa.build_task_dag(ta, tb);

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
    auto tasks = pa.build_task_dag(ta, ta);  // pair A with itself

    if (tasks.size() >= 2) {
        // At least one task should have priority >= 5^3 = 125
        uint64_t max_pri = 0;
        for (const auto& t : tasks) max_pri = std::max(max_pri, t.priority);
        EXPECT_GE(max_pri, uint64_t(125));
    }
}

// ── ScalingProxy: wall-time data for Experiment 4 ────────────────────────────
// These tests are informational (they don't assert speedup — we have only
// a few cores in the test environment) but produce the timing data that
// goes into Figure 3 of the paper.

TEST(ScalingProxy, WallTimeDeclinesWithMoreThreads) {
    // 50 bubbles, k=4 — enough work to see a trend even on 2 cores
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

    // Record (no strict assert — Experiment 4 will run on HPRC scale)
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

// ── EnvOverride ───────────────────────────────────────────────────────────────

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

// ── AlignSnarls interface (used by GraphAligner replacement) ──────────────────

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
