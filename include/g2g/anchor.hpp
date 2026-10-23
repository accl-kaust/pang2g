#pragma once
/**
 * include/g2g/anchor.hpp
 * ======================
 * Anchor node detection between two variation graphs.
 *
 * Phase 2 finding: _find_anchors accounts for 39% of runtime at n=500
 * via O(|VA|·|VB|) naïve pairwise scan.  At HPRC scale |V|~10^7,
 * this is 10^14 comparisons - completely infeasible.
 *
 * This module replaces the naïve scan with a minimizer k-mer index
 * (sliding-window sketching, as in minimap2):
 *
 *   1. Build a minimizer index over graph A's backbone node sequences.
 *   2. Query each backbone node in graph B: candidate anchors are pairs
 *      sharing ≥1 minimizer.
 *   3. Refine candidates with full seq_similarity; keep pairs ≥ min_sim.
 *
 * The index uses 32-bit hash keys so the lookup fits in L2 cache for
 * typical backbone sizes (~10^5–10^6 nodes).
 *
 * Phase 3 task estimate: 2 weeks (§3.3, item 5).
 */

#include <g2g/types.hpp>
#include <g2g/graph.hpp>

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace g2g {

struct AnchorParams {
    float    min_sim    = 0.90f;  // minimum similarity to call a pair an anchor
    uint8_t  kmer_len   = 15;     // k-mer length for minimizer sketching
    uint8_t  window     = 10;     // minimizer window size (as in minimap2)
};

/**
 * MinimizerIndex
 * Compact k-mer sketch over one graph's backbone node sequences.
 * Backbone = nodes with in-degree ≤ 1 AND out-degree ≤ 1.
 */
class MinimizerIndex {
public:
    explicit MinimizerIndex(const VariationGraph& vg,
                            const AnchorParams& params = {});

    /**
     * For a query sequence, return the set of node IDs in the indexed
     * graph that share ≥1 minimizer with the query.
     */
    std::vector<NodeId> query(std::string_view query_seq) const;

    uint32_t indexed_nodes() const { return static_cast<uint32_t>(backbone_.size()); }

private:
    const VariationGraph& vg_;
    AnchorParams          params_;

    std::vector<NodeId>   backbone_;   // backbone node ids

    // minimizer hash → list of node ids that produce it
    std::unordered_map<uint32_t, std::vector<NodeId>> index_;

    void build();
    std::vector<uint32_t> minimizers(std::string_view seq) const;
    static uint32_t hash_kmer(std::string_view kmer);
};

/**
 * Primary entry point.
 * Returns anchor pairs (na, nb) with sim ≥ params.min_sim,
 * deduplicated (each node used at most once, greedy best-first).
 *
 * Uses MinimizerIndex when |backbone_a| · |backbone_b| > 10^8,
 * falls back to naïve pairwise scan for small graphs (unit tests).
 */
std::vector<AnchorMatch> find_anchors(const VariationGraph& ga,
                                       const VariationGraph& gb,
                                       const AnchorParams& params = {});

} // namespace g2g
