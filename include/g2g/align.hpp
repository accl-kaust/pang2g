#pragma once
/**
 * include/g2g/align.hpp
 * =====================
 * Graph-to-graph alignment — sequential driver and snarl pairing.
 *
 * Pipeline (GraphAligner::align, mirrored by ParallelAligner::align):
 *   1. find_anchors()        — backbone node matching (anchor.hpp)
 *   2. SnarlDecomposer       — snarl tree construction (decompose.hpp)
 *   3. pair_snarls()         — snarl pairing (boundary-key or anchor-interval)
 *   4. score_snarl_pair()    — W_seq + W_topo + W_struct (LAPJV allele matching)
 *   5. build_deltas()        — delta graphs: unmatched alleles of paired snarls
 *                              plus every allele of an unpaired snarl
 *
 * The parallel runtime (parallel.hpp) replaces the per-pair loop of step 4 with
 * a priority-dispatch task DAG; steps 1-3 and 5 are shared.
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

/**
 * How snarls of A are paired with snarls of B before the per-pair matching.
 *
 *  Boundary  — identical (source, sink) boundary sequences; repeated keys are
 *              consumed positionally (FIFO).  O(|TA|+|TB|).  Default; this is
 *              the rule all reported scaling/ablation numbers were produced with.
 *  Interval  — anchor-guided: snarls are grouped into the intervals between
 *              consecutive shared anchors and paired within each interval by
 *              an order-preserving dynamic program.  Contains pairing errors to
 *              a single interval.  Near-linear when anchors are dense.
 */
enum class PairingMode { Boundary, Interval };

struct AlignerParams {
    float alpha         = 1.0f;   // W_seq weight
    float beta          = 0.5f;   // W_topo weight
    float gamma         = 2.0f;   // W_struct weight
    float gap           = 0.5f;   // gap penalty per unmatched allele

    PairingMode  pairing = PairingMode::Boundary;

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
    double dp_seconds        = 0.0;  // snarl pairing + per-pair DP (align_snarls)
};

class GraphAligner {
public:
    explicit GraphAligner(const AlignerParams& params = {})
        : params_(params) {}

    /** Run the full graph-to-graph alignment, single-threaded. */
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
};

/**
 * Boundary-key snarl correspondence (PairingMode::Boundary).
 *
 * Two snarls correspond when their boundary (source, sink) sequences are
 * byte-identical.  Pairing is injective and resolved in snarl order, so
 * repeated boundary sequences fall back to a positional (FIFO) match.  Runs in
 * O(|TA| + |TB|).  A snarl present in only one graph whose key repeats shifts
 * every later snarl with that key by one position — the cascade that
 * PairingMode::Interval prevents.
 */
std::vector<std::pair<uint32_t, uint32_t>> correspond_snarls(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb);

/**
 * Anchor-interval snarl correspondence (PairingMode::Interval).
 *
 * The anchors are chained in topological order (longest chain increasing in
 * both graphs); the snarls between two consecutive anchors are paired by an
 * order-preserving dynamic program whose match score is the fraction of
 * identical allele sequences plus the boundary-node similarity.  Snarls with
 * no identical allele are never paired (their alleles all go to the delta).
 * Cost: sum over intervals of |A_i|·|B_i|; an interval exceeding 2^24 cells
 * falls back to the boundary rule.  With no anchors the whole graph is one
 * interval.
 */
std::vector<std::pair<uint32_t, uint32_t>> correspond_snarls_interval(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb,
    const std::vector<AnchorMatch>& anchors);

/** Dispatch on `mode`; shared by the sequential and parallel aligners. */
std::vector<std::pair<uint32_t, uint32_t>> pair_snarls(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb,
    const std::vector<AnchorMatch>& anchors,
    PairingMode mode);

/**
 * Fill result.a_only_alleles / b_only_alleles (and unpaired_snarls_a/b) from
 * result.snarl_alignments: the alleles LAPJV left unmatched inside paired
 * snarls, plus every allele of a snarl that was not paired at all.  Shared by
 * GraphAligner and ParallelAligner so both emit byte-identical deltas.
 */
void build_deltas(const VariationGraph& ga,
                  const VariationGraph& gb,
                  const SnarlTree& ta,
                  const SnarlTree& tb,
                  AlignmentResult& result);

} // namespace g2g
