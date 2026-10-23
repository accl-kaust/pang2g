#pragma once
/**
 * include/g2g/sequence.hpp
 * ========================
 * String-level computation - C++ port of sequence.py.
 *
 * Phase 3 optimisation targets (per Phase 2 cProfile findings):
 *   - edit_distance  accounts for 70% of total runtime via 2.1M calls.
 *     The inner Wagner-Fischer DP loop is a prime candidate for
 *     SIMD vectorisation (SSE4.2 / AVX2 PCMPEQB / VPCMPEQB).
 *   - seq_similarity is a thin wrapper; keep it inline.
 *   - path_sequence should concatenate into a pre-allocated buffer.
 */

#include <g2g/graph.hpp>
#include <g2g/decompose.hpp>

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>
#include <g2g/span.hpp>   // portable g2g::span (falls back on GCC 9)

namespace g2g {

/**
 * Wagner-Fischer edit distance.
 * Space O(n), Time O(m*n).
 * TODO Phase 3: add SIMD fast-path for short strings (≤16 bytes → SSE4.2).
 */
int edit_distance(std::string_view a, std::string_view b);

/**
 * Normalised similarity in [0, 1]: 1 - edit_dist / max(|a|, |b|).
 * Both empty → 1.0.
 */
inline float seq_similarity(std::string_view a, std::string_view b) {
    if (a.empty() && b.empty()) return 1.0f;
    size_t denom = std::max(a.size(), b.size());
    return 1.0f - static_cast<float>(edit_distance(a, b)) / static_cast<float>(denom);
}

/**
 * Concatenate node sequences along a path into out_buf.
 * Appends to out_buf (does not clear first) so callers can reuse allocations.
 */
void path_sequence(const VariationGraph& vg,
                   g2g::span<const NodeId> path,
                   std::string& out_buf);

/** Convenience overload returning a new string. */
inline std::string path_sequence(const VariationGraph& vg,
                                  g2g::span<const NodeId> path) {
    std::string buf;
    path_sequence(vg, path, buf);
    return buf;
}

/**
 * Compute allele sequences for a snarl: inner nodes only
 * (source and sink are excluded as shared boundary nodes).
 * Fills out - one string per allele - reusing buffer capacity.
 * Requires the SnarlTree that owns the allele-path store for this snarl.
 */
void allele_sequences(const VariationGraph& vg,
                      const Snarl& s,
                      const SnarlTree& tree,
                      std::vector<std::string>& out);

} // namespace g2g
