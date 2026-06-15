#pragma once
/**
 * include/g2g/parallel.hpp
 * ========================
 * Phase 4 — Work-stealing parallel runtime for the snarl-pair DP.
 *
 * OpenMP implementation (replaces oneTBB).
 *
 * Design
 * ------
 * The snarl-pair alignment has a natural task DAG: each pair (SA, SB) depends
 * on all of its child pairs (bottom-up DP).  Sibling pairs are fully
 * independent and can run in parallel.
 *
 * We implement the scheduler with:
 *   1. SnarlPairTask  — A unit of work: score one (SA_idx, SB_idx) pair.
 *   2. A std::priority_queue protected by a mutex as the ready queue,
 *      ordered by subtree_work descending (steal-largest).
 *   3. OpenMP threads pulling tasks from the ready queue.
 *
 * Thread count
 * ------------
 * Defaults to omp_get_max_threads().  Set G2G_NUM_THREADS in the
 * environment or pass explicitly to ParallelAligner to override.
 */

#include <g2g/types.hpp>
#include <g2g/decompose.hpp>
#include <g2g/align.hpp>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace g2g {

// ── Task record ───────────────────────────────────────────────────────────────

struct SnarlPairTask {
    uint32_t ia;           // index into SnarlTree A
    uint32_t ib;           // index into SnarlTree B
    uint64_t priority;     // subtree_work — scheduler hint

    int32_t  parent_idx;   // -1 if root
    std::atomic<int32_t> dep_count{0};

    std::unique_ptr<SnarlAlignment> result;

    SnarlPairTask() = default;
    SnarlPairTask(const SnarlPairTask&) = delete;
    SnarlPairTask& operator=(const SnarlPairTask&) = delete;
    SnarlPairTask(SnarlPairTask&& o) noexcept
        : ia(o.ia), ib(o.ib), priority(o.priority)
        , parent_idx(o.parent_idx)
        , dep_count(o.dep_count.load(std::memory_order_relaxed))
        , result(std::move(o.result))
    {}
    SnarlPairTask& operator=(SnarlPairTask&&) = delete;
};

// ── Scheduler statistics ──────────────────────────────────────────────────────

struct RunStats {
    uint64_t tasks_executed{0};
    uint64_t tasks_ready_peak{0};
    double   wall_seconds{0.0};
    uint32_t num_threads{0};

    double   parallel_efficiency{0.0};
    double   sequential_seconds{0.0};
};

// ── Parallel aligner ──────────────────────────────────────────────────────────

class ParallelAligner {
public:
    explicit ParallelAligner(const AlignerParams& params = {},
                              uint32_t num_threads = 0);

    struct Result {
        std::vector<SnarlAlignment> alignments;
        RunStats                    stats;
    };

    Result align_snarls(const VariationGraph& ga,
                        const VariationGraph& gb,
                        const SnarlTree& ta,
                        const SnarlTree& tb,
                        const std::vector<AnchorMatch>& anchors) const;

    AlignmentResult align(const VariationGraph& ga,
                          const VariationGraph& gb) const;

    uint32_t    num_threads() const { return num_threads_; }
    RunStats    last_stats()  const { return last_stats_; }

    std::vector<SnarlPairTask> build_task_dag(
        const SnarlTree& ta,
        const SnarlTree& tb) const;

private:
    AlignerParams    params_;
    uint32_t         num_threads_;
    mutable RunStats last_stats_;

    SnarlAlignment execute_task(
        const VariationGraph& ga,
        const VariationGraph& gb,
        const SnarlTree& ta,
        const SnarlTree& tb,
        const SnarlPairTask& task) const;
};

} // namespace g2g
