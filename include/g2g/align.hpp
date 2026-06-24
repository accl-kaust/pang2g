#pragma once
/**
 * include/g2g/align.hpp
 * =====================
 * Graph-to-graph alignment — C++ port of GraphAligner from align.py.
 *
 * Algorithm is identical to the Python prototype:
 *   1. find_anchors()       — backbone node matching (anchor.hpp)
 *   2. SnarlDecomposer      — snarl tree construction (decompose.hpp)
 *   3. _align_snarls()      — greedy snarl pairing + per-pair DP
 *   4. _score_snarl_pair()  — W_seq + W_topo + W_struct
 *   5. _build_deltas()      — fill AlignmentResult delta fields
 *
 * Phase 4 hook: GraphAligner::align() is the entry point that Phase 4
 * will parallelise by replacing the sequential snarl loop with a
 * work-stealing task DAG (see align.cpp for the TODO marker).
 */

#include <g2g/types.hpp>
#include <g2g/graph.hpp>
#include <g2g/decompose.hpp>
#include <g2g/anchor.hpp>
#include <g2g/matching.hpp>

#include <cstdint>
#include <utility>
#include <vector>

namespace g2g {

struct AlignerParams {
    float alpha         = 1.0f;   // W_seq weight
    float beta          = 0.5f;   // W_topo weight
    float gamma         = 2.0f;   // W_struct weight
    float gap           = 0.5f;   // gap penalty per unmatched allele

    AnchorParams anchor;
    MatchParams  match;
};

// ── Per-phase wall-clock breakdown ────────────────────────────────────────────
// Populated by the aligners so the CLI can emit a `[timing]` line under
// --verbose-timing.  (I/O time is measured separately in main, since GFA
// reading happens before the aligner is invoked.)
struct PhaseTimes {
    double index_seconds     = 0.0;  // anchor detection / minimizer index
    double decompose_seconds = 0.0;  // snarl-tree decomposition (both graphs)
    double dp_seconds        = 0.0;  // snarl-pair DP (align_snarls)
};

class GraphAligner {
public:
    explicit GraphAligner(const AlignerParams& params = {})
        : params_(params) {}

    /**
     * Run full graph-to-graph alignment.
     * Single-threaded in Phase 3; Phase 4 replaces inner loop with
     * work-stealing scheduler.
     */
    AlignmentResult align(const VariationGraph& ga,
                          const VariationGraph& gb) const;

    // Per-phase wall-clock breakdown of the most recent align() call.
    PhaseTimes last_phase_times() const { return phase_times_; }

private:
    AlignerParams params_;
    mutable PhaseTimes phase_times_;

    std::vector<SnarlAlignment> align_snarls(
        const VariationGraph& ga,
        const VariationGraph& gb,
        const SnarlTree& ta,
        const SnarlTree& tb,
        const std::vector<AnchorMatch>& anchors) const;

    SnarlAlignment score_snarl_pair(
        const VariationGraph& ga,
        const VariationGraph& gb,
        const SnarlTree& ta,
        const SnarlTree& tb,
        uint32_t ia,
        uint32_t ib) const;

    void build_deltas(const VariationGraph& ga,
                      const VariationGraph& gb,
                      AlignmentResult& result,
                      const SnarlTree& ta,
                      const SnarlTree& tb) const;
};

/**
 * Establish the snarl-pair correspondence between two graphs.
 *
 * Two snarls correspond when their boundary (source, sink) sequences match —
 * the same signal the per-pair score maximises (W_seq) and the structure the
 * pangenome backbone provides (shared anchor nodes between A and B).  Pairing
 * is injective and resolved in snarl order, so repeated boundary sequences fall
 * back to a stable positional match.  Runs in O(|TA| + |TB|).
 *
 * The parallel scheduler pairs through this function (replacing an earlier
 * heuristic that paired by subtree work — i.e. snarl size — and mismatched
 * snarls whenever their sizes differed, yielding alignments inconsistent with
 * the sequential aligner).  On the structured graphs the project's generators
 * emit, this reproduces the same correspondence the sequential aligner's
 * best-score search finds, so the two runtimes agree — while staying O(n)
 * rather than the sequential O(|TA|·|TB|) scan.
 */
std::vector<std::pair<uint32_t, uint32_t>> correspond_snarls(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb);

} // namespace g2g
