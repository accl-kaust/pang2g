/**
 * src/parallel.cpp
 * ================
 * Phase 4 — Work-stealing parallel runtime (OpenMP implementation).
 *
 * Replaces oneTBB with:
 *   - OpenMP threads (#pragma omp parallel)
 *   - std::priority_queue + std::mutex as the steal-largest ready queue
 *   - std::condition_variable for idle-wait instead of spin-yield
 *
 * The steal-largest semantics are preserved: workers pop the highest-
 * priority (largest subtree_work) task from the shared ready queue.
 * OpenMP's own work-stealing applies to the outer parallel region;
 * our explicit priority queue ensures large tasks are picked first
 * regardless of which thread happens to wake up.
 *
 * G2G_STEAL_RANDOM ablation: when compiled with -DG2G_STEAL_RANDOM=1,
 * the priority queue is replaced with a FIFO std::queue (tasks popped
 * in insertion order rather than by priority).
 *
 * G2G_STATIC_PARTITION ablation: when compiled with -DG2G_STATIC_PARTITION=1,
 * the shared queue and stealing are removed entirely — each thread owns a
 * fixed (cyclic) slice of the tasks and runs only those.  This exposes the
 * per-snarl load imbalance that work-stealing hides; results are unchanged.
 */

#include <g2g/parallel.hpp>
#include <g2g/sequence.hpp>
#include <g2g/matching.hpp>

#ifdef G2G_HAVE_OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cassert>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <mutex>
#include <numeric>
#include <queue>
#include <stdexcept>

namespace g2g {

// ── Helpers ───────────────────────────────────────────────────────────────────

static uint32_t detect_threads() {
    if (const char* env = std::getenv("G2G_NUM_THREADS")) {
        int n = std::atoi(env);
        if (n > 0) return static_cast<uint32_t>(n);
    }
#ifdef G2G_HAVE_OPENMP
    return static_cast<uint32_t>(omp_get_max_threads());
#else
    uint32_t hw = std::thread::hardware_concurrency();
    return hw ? hw : 1u;
#endif
}

// ── ParallelAligner constructor ───────────────────────────────────────────────

ParallelAligner::ParallelAligner(const AlignerParams& params, uint32_t num_threads)
    : params_(params)
    , num_threads_(num_threads == 0 ? detect_threads() : num_threads)
{}

// ── Task execution (identical to sequential score_snarl_pair) ─────────────────

SnarlAlignment ParallelAligner::execute_task(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb,
    const SnarlPairTask& task) const
{
    const Snarl& sa = ta.snarls[task.ia];
    const Snarl& sb = tb.snarls[task.ib];

    SnarlAlignment aln{};
    aln.snarl_a_idx = task.ia;
    aln.snarl_b_idx = task.ib;

    // W_seq
    float src_sim  = seq_similarity(ga.seq(sa.source), gb.seq(sb.source));
    float sink_sim = seq_similarity(ga.seq(sa.sink),   gb.seq(sb.sink));
    aln.w_seq = params_.alpha * (src_sim + sink_sim) / 2.f;

    // W_topo
    aln.w_topo = params_.beta * 1.f;

    // W_struct
    std::vector<std::string> seqs_a_str, seqs_b_str;
    allele_sequences(ga, sa, ta, seqs_a_str);
    allele_sequences(gb, sb, tb, seqs_b_str);

    std::vector<std::string_view> seqs_a_sv, seqs_b_sv;
    seqs_a_sv.reserve(seqs_a_str.size());
    seqs_b_sv.reserve(seqs_b_str.size());
    for (const auto& s : seqs_a_str) seqs_a_sv.emplace_back(s);
    for (const auto& s : seqs_b_str) seqs_b_sv.emplace_back(s);

    MatchParams mp;
    mp.gap_cost = params_.gap;
    mp.min_sim  = params_.match.min_sim;

    MatchResult mr = allele_match(seqs_a_sv, seqs_b_sv, mp);
    aln.w_struct       = params_.gamma * mr.score;
    aln.allele_matches = mr.matches;
    aln.a_only         = mr.a_unmatched;
    aln.b_only         = mr.b_unmatched;

    aln.score = aln.w_seq + aln.w_topo + aln.w_struct;
    return aln;
}

// ── Task DAG construction ─────────────────────────────────────────────────────

std::vector<SnarlPairTask> ParallelAligner::build_task_dag(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb) const
{
    const uint32_t na = static_cast<uint32_t>(ta.snarls.size());

    std::vector<int32_t> task_idx_for_ia(na, -1);

    // Pair snarls by boundary correspondence (shared with the sequential
    // aligner's intent) rather than by subtree size.  subtree_work is used
    // only as the scheduler priority below.
    std::vector<std::pair<uint32_t, uint32_t>> corr =
        correspond_snarls(ga, gb, ta, tb);

    std::vector<SnarlPairTask> tasks;
    tasks.reserve(corr.size());

    for (const auto& [ia, ib] : corr) {
        uint32_t task_idx = static_cast<uint32_t>(tasks.size());
        task_idx_for_ia[ia] = static_cast<int32_t>(task_idx);

        SnarlPairTask& t = tasks.emplace_back();
        t.ia         = ia;
        t.ib         = ib;
        t.priority   = ta.snarls[ia].subtree_work + tb.snarls[ib].subtree_work;
        t.parent_idx = -1;
        t.dep_count.store(0, std::memory_order_relaxed);
    }

    // Wire parent/child dependencies from snarl-tree nesting
    std::vector<uint32_t> parent_of_ia(na, UINT32_MAX);
    for (uint32_t i = 0; i < na; ++i) {
        for (uint32_t c = ta.snarls[i].children_begin;
             c < ta.snarls[i].children_end; ++c) {
            if (c < na) parent_of_ia[c] = i;
        }
    }

    for (uint32_t ti = 0; ti < static_cast<uint32_t>(tasks.size()); ++ti) {
        uint32_t ia  = tasks[ti].ia;
        uint32_t par = parent_of_ia[ia];
        if (par != UINT32_MAX && task_idx_for_ia[par] >= 0) {
            uint32_t par_ti = static_cast<uint32_t>(task_idx_for_ia[par]);
            tasks[ti].parent_idx = static_cast<int32_t>(par_ti);
            tasks[par_ti].dep_count.fetch_add(1, std::memory_order_relaxed);
        }
    }

    return tasks;
}

// ── Ready-queue entry ─────────────────────────────────────────────────────────

struct ReadyEntry {
    uint32_t task_idx;
    uint64_t priority;
};

struct CmpPriority {
    bool operator()(const ReadyEntry& a, const ReadyEntry& b) const {
        return a.priority < b.priority;   // max-heap
    }
};

// ── Main parallel scheduler (OpenMP) ──────────────────────────────────────────

ParallelAligner::Result ParallelAligner::align_snarls(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb,
    const std::vector<AnchorMatch>& /*anchors*/) const
{
    using Clock = std::chrono::steady_clock;
    auto t0 = Clock::now();

    std::vector<SnarlPairTask> tasks = build_task_dag(ga, gb, ta, tb);
    const uint32_t ntasks = static_cast<uint32_t>(tasks.size());

    // Per-task result slots
    std::vector<SnarlAlignment>    results(ntasks);
    std::vector<std::atomic<bool>> result_ready(ntasks);
    for (auto& f : result_ready) f.store(false, std::memory_order_relaxed);

    std::atomic<uint64_t> tasks_executed{0};
    std::atomic<uint64_t> steals{0};

    // Per-thread busy time (sum of execute_task durations).  Each thread writes
    // only its own slot, so no synchronisation is needed; we reduce after the
    // parallel region to derive the load-imbalance factor.
    std::vector<double> thread_busy(num_threads_, 0.0);
    uint64_t ready_peak = 0;   // guarded by queue_mtx (and the pre-launch seed)

#ifdef G2G_STATIC_PARTITION
    // ── Static work-partition scheduler (ablation) ────────────────────────────
    // No shared queue and no stealing: thread `tid` of a P-thread team owns
    // exactly the tasks with (ti % P == tid) and runs only those, re-scanning its
    // own set until each owned task's dependencies are met.  Dependencies are
    // still honoured, so the results are byte-identical to the work-stealing path
    // — but a thread that draws several heavy snarls cannot shed them, so the DP
    // phase inherits the per-snarl load imbalance that work-stealing removes.  By
    // construction no task ever crosses its static owner, so `steals` stays 0 and
    // `imbalance` reports the static skew.
    ready_peak = ntasks;                 // every task is "owned" from the start
    auto worker = [&]() {
        const uint32_t P =
#ifdef G2G_HAVE_OPENMP
            static_cast<uint32_t>(omp_get_num_threads());
#else
            1u;
#endif
        const uint32_t tid =
#ifdef G2G_HAVE_OPENMP
            static_cast<uint32_t>(omp_get_thread_num());
#else
            0u;
#endif
        if (P == 0 || tid >= num_threads_) return;
        std::vector<uint32_t> mine;
        for (uint32_t ti = tid; ti < ntasks; ti += P) mine.push_back(ti);
        std::vector<char> finished(mine.size(), 0);
        uint32_t remaining = static_cast<uint32_t>(mine.size());
        while (remaining > 0) {
            bool progressed = false;
            for (std::size_t k = 0; k < mine.size(); ++k) {
                if (finished[k]) continue;
                const uint32_t ti = mine[k];
                if (tasks[ti].dep_count.load(std::memory_order_acquire) != 0)
                    continue;            // a child (maybe on another thread) is pending
                auto exec_t0 = Clock::now();
                results[ti] = execute_task(ga, gb, ta, tb, tasks[ti]);
                thread_busy[tid] +=
                    std::chrono::duration<double>(Clock::now() - exec_t0).count();
                result_ready[ti].store(true, std::memory_order_release);
                tasks_executed.fetch_add(1, std::memory_order_relaxed);
                finished[k] = 1; --remaining; progressed = true;
                if (tasks[ti].parent_idx >= 0) {
                    uint32_t pi = static_cast<uint32_t>(tasks[ti].parent_idx);
                    tasks[pi].dep_count.fetch_sub(1, std::memory_order_acq_rel);
                }
            }
            if (!progressed && remaining > 0)
                std::this_thread::yield();          // wait on another thread's child
        }
    };
#else
    std::atomic<uint32_t> done_count{0};
    // Shared ready queue protected by mutex + condvar
#ifdef G2G_STEAL_RANDOM
    std::queue<ReadyEntry> ready_queue;
#else
    std::priority_queue<ReadyEntry, std::vector<ReadyEntry>, CmpPriority> ready_queue;
#endif
    std::mutex              queue_mtx;
    std::condition_variable queue_cv;

    // Seed with all immediately-ready tasks
    for (uint32_t ti = 0; ti < ntasks; ++ti) {
        if (tasks[ti].dep_count.load(std::memory_order_relaxed) == 0) {
            ready_queue.push({ti, tasks[ti].priority});
        }
    }
    ready_peak = ready_queue.size();

    // Worker function: pull tasks from the ready queue until all done
    auto worker = [&]() {
        const uint32_t tid =
#ifdef G2G_HAVE_OPENMP
            static_cast<uint32_t>(omp_get_thread_num());
#else
            0u;
#endif
        while (true) {
            ReadyEntry entry{};
            bool got_task = false;

            {
                std::unique_lock<std::mutex> lk(queue_mtx);
                // Wait until there's a task or all tasks are done
                queue_cv.wait(lk, [&]() {
                    return !ready_queue.empty()
                        || done_count.load(std::memory_order_acquire) >= ntasks;
                });

                if (done_count.load(std::memory_order_acquire) >= ntasks)
                    return;

                if (!ready_queue.empty()) {
#ifdef G2G_STEAL_RANDOM
                    entry = ready_queue.front();
                    ready_queue.pop();
#else
                    entry = ready_queue.top();
                    ready_queue.pop();
#endif
                    got_task = true;
                }
            }

            if (!got_task) continue;

            uint32_t ti = entry.task_idx;
            const SnarlPairTask& task = tasks[ti];

            // Steal accounting: a task whose static round-robin owner
            // (ti % P) differs from the executing thread was redistributed
            // by the dynamic scheduler rather than handled by a static
            // partition.  (For the current flat snarl DAGs this captures
            // initial work distribution; it also covers dynamically-enqueued
            // child tasks once nested snarls are produced.)
            if (static_cast<uint32_t>(ti % num_threads_) != tid)
                steals.fetch_add(1, std::memory_order_relaxed);

            // Execute (timed, to derive per-thread load imbalance)
            auto exec_t0 = Clock::now();
            results[ti] = execute_task(ga, gb, ta, tb, task);
            thread_busy[tid] +=
                std::chrono::duration<double>(Clock::now() - exec_t0).count();
            result_ready[ti].store(true, std::memory_order_release);

            tasks_executed.fetch_add(1, std::memory_order_relaxed);
            uint32_t new_done = done_count.fetch_add(1, std::memory_order_release) + 1;

            // Notify parent
            if (task.parent_idx >= 0) {
                uint32_t pi = static_cast<uint32_t>(task.parent_idx);
                int32_t remaining = tasks[pi].dep_count.fetch_sub(
                    1, std::memory_order_acq_rel) - 1;
                if (remaining == 0) {
                    {
                        std::lock_guard<std::mutex> lk(queue_mtx);
                        ready_queue.push({pi, tasks[pi].priority});
                        ready_peak = std::max(
                            ready_peak,
                            static_cast<uint64_t>(ready_queue.size()));
                    }
                    queue_cv.notify_one();
                }
            }

            // If all done, wake everyone so they can exit.  done_count is
            // updated outside queue_mtx (line above), so the notify must be
            // issued while holding the lock: otherwise a worker that read the
            // pre-increment done_count in its wait predicate could go to sleep
            // *after* this notify fires and never be woken again (lost wakeup →
            // the OpenMP region never joins and the run hangs).
            if (new_done >= ntasks) {
                std::lock_guard<std::mutex> lk(queue_mtx);
                queue_cv.notify_all();
            }
        }
    };
#endif // G2G_STATIC_PARTITION scheduler select

    // Launch OpenMP threads
#ifdef G2G_HAVE_OPENMP
    #pragma omp parallel num_threads(num_threads_)
    {
        worker();
    }
#else
    // Fallback: single-threaded
    worker();
#endif

    auto t1 = Clock::now();
    double wall = std::chrono::duration<double>(t1 - t0).count();

    std::vector<SnarlAlignment> alignments;
    alignments.reserve(ntasks);
    for (uint32_t ti = 0; ti < ntasks; ++ti)
        alignments.push_back(std::move(results[ti]));

    // Load-imbalance factor: max thread busy-time over the mean.  Equals 1.0
    // for a perfectly balanced run and approaches P in the worst case.
    double busy_total = std::accumulate(thread_busy.begin(),
                                        thread_busy.end(), 0.0);
    double busy_max   = thread_busy.empty()
        ? 0.0
        : *std::max_element(thread_busy.begin(), thread_busy.end());
    double busy_mean  = num_threads_ > 0
        ? busy_total / static_cast<double>(num_threads_)
        : 0.0;

    RunStats stats;
    stats.tasks_executed   = tasks_executed.load();
    stats.tasks_ready_peak = ready_peak;
    stats.steals           = steals.load();
    stats.imbalance        = busy_mean > 0.0 ? busy_max / busy_mean : 1.0;
    stats.wall_seconds     = wall;
    stats.num_threads      = num_threads_;

    last_stats_ = stats;
    return {std::move(alignments), stats};
}

// ── Full pipeline ─────────────────────────────────────────────────────────────

AlignmentResult ParallelAligner::align(
    const VariationGraph& ga,
    const VariationGraph& gb) const
{
    AlignmentResult result;

    using Clock = std::chrono::steady_clock;
    auto secs = [](Clock::time_point a, Clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    };

    auto t_idx0 = Clock::now();
    result.anchors = find_anchors(ga, gb, params_.anchor);
    auto t_idx1 = Clock::now();

    SnarlTree ta = SnarlDecomposer(ga).decompose();
    SnarlTree tb = SnarlDecomposer(gb).decompose();
    auto t_dec1 = Clock::now();

    auto [alignments, stats] = align_snarls(ga, gb, ta, tb, result.anchors);
    result.snarl_alignments  = std::move(alignments);

    stats.phases.index_seconds     = secs(t_idx0, t_idx1);
    stats.phases.decompose_seconds = secs(t_idx1, t_dec1);
    stats.phases.dp_seconds        = stats.wall_seconds;  // align_snarls wall
    last_stats_              = stats;

    for (const auto& a : result.anchors)
        result.total_score += params_.alpha * a.sim;
    for (const auto& sa : result.snarl_alignments)
        result.total_score += sa.score;

    // Delta graph construction
    auto node_name = [](const VariationGraph& vg, NodeId n) -> std::string {
        if (n < vg.node_names.size()) return vg.node_names[n];
        return std::to_string(n);
    };
    auto make_snarl_id = [&](const VariationGraph& vg, const Snarl& s) {
        return node_name(vg, s.source) + ">" + node_name(vg, s.sink);
    };
    auto inner_seq = [](const VariationGraph& vg,
                         g2g::span<const NodeId> full_path) -> std::string {
        std::string seq;
        if (full_path.size() > 2) {
            for (NodeId n : full_path.subspan(1, full_path.size() - 2)) {
                auto sv = vg.seq(n);
                seq.append(sv.data(), sv.size());
            }
        }
        return seq;
    };

    for (const auto& aln : result.snarl_alignments) {
        if (aln.snarl_a_idx >= ta.snarls.size()) continue;
        if (aln.snarl_b_idx >= tb.snarls.size()) continue;
        const Snarl& sa = ta.snarls[aln.snarl_a_idx];
        const Snarl& sb = tb.snarls[aln.snarl_b_idx];
        for (uint16_t idx : aln.a_only) {
            auto path = ta.allele_path(sa.alleles_begin, idx);
            DeltaAllele da;
            da.snarl_id = make_snarl_id(ga, sa);
            da.path     = std::vector<NodeId>(path.begin(), path.end());
            da.sequence = inner_seq(ga, path);
            result.a_only_alleles.push_back(std::move(da));
        }
        for (uint16_t idx : aln.b_only) {
            auto path = tb.allele_path(sb.alleles_begin, idx);
            DeltaAllele db;
            db.snarl_id = make_snarl_id(gb, sb);
            db.path     = std::vector<NodeId>(path.begin(), path.end());
            db.sequence = inner_seq(gb, path);
            result.b_only_alleles.push_back(std::move(db));
        }
    }

    return result;
}

} // namespace g2g
