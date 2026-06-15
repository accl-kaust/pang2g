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
    const SnarlTree& ta,
    const SnarlTree& tb) const
{
    const uint32_t na = static_cast<uint32_t>(ta.snarls.size());
    const uint32_t nb = static_cast<uint32_t>(tb.snarls.size());

    std::vector<bool>    used_b(nb, false);
    std::vector<int32_t> task_idx_for_ia(na, -1);

    std::vector<SnarlPairTask> tasks;
    tasks.reserve(std::min(na, nb));

    for (uint32_t ia = 0; ia < na; ++ia) {
        int32_t  best_ib  = -1;
        uint64_t best_pri = 0;
        for (uint32_t ib = 0; ib < nb; ++ib) {
            if (used_b[ib]) continue;
            uint64_t pri = ta.snarls[ia].subtree_work + tb.snarls[ib].subtree_work;
            if (best_ib < 0 || pri > best_pri) {
                best_pri = pri;
                best_ib  = static_cast<int32_t>(ib);
            }
        }
        if (best_ib < 0) break;

        used_b[static_cast<uint32_t>(best_ib)] = true;

        uint32_t task_idx = static_cast<uint32_t>(tasks.size());
        task_idx_for_ia[ia] = static_cast<int32_t>(task_idx);

        SnarlPairTask& t = tasks.emplace_back();
        t.ia         = ia;
        t.ib         = static_cast<uint32_t>(best_ib);
        t.priority   = best_pri;
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

    std::vector<SnarlPairTask> tasks = build_task_dag(ta, tb);
    const uint32_t ntasks = static_cast<uint32_t>(tasks.size());

    // Per-task result slots
    std::vector<SnarlAlignment>    results(ntasks);
    std::vector<std::atomic<bool>> result_ready(ntasks);
    for (auto& f : result_ready) f.store(false, std::memory_order_relaxed);

    std::atomic<uint64_t> tasks_executed{0};
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

    // Worker function: pull tasks from the ready queue until all done
    auto worker = [&]() {
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

            // Execute
            results[ti] = execute_task(ga, gb, ta, tb, task);
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
                    }
                    queue_cv.notify_one();
                }
            }

            // If all done, wake everyone so they can exit
            if (new_done >= ntasks) {
                queue_cv.notify_all();
            }
        }
    };

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

    RunStats stats;
    stats.tasks_executed   = tasks_executed.load();
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

    result.anchors = find_anchors(ga, gb, params_.anchor);

    SnarlTree ta = SnarlDecomposer(ga).decompose();
    SnarlTree tb = SnarlDecomposer(gb).decompose();

    auto [alignments, stats] = align_snarls(ga, gb, ta, tb, result.anchors);
    result.snarl_alignments  = std::move(alignments);
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
                         std::span<const NodeId> full_path) -> std::string {
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
