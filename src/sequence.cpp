/**
 * src/sequence.cpp
 * ================
 * C++ port of sequence.py.
 *
 * edit_distance: Wagner-Fischer DP, space O(n).
 * TODO Phase 3 SIMD: for sequences ≤ 16 bytes (most inner-node seqs),
 *   use SSE4.2 PCMPEQB to compare 16 chars in one instruction.
 *   Speedup estimate: 8–16× on x86 (see §2.5.1 finding 1 in the paper).
 *   Gate behind #ifdef __SSE4_2__ so the scalar path remains the fallback.
 */

#include <g2g/sequence.hpp>

#include <algorithm>
#include <cassert>
#include <vector>

namespace g2g {

int edit_distance(std::string_view a, std::string_view b) {
    const size_t m = a.size(), n = b.size();
    if (m == 0) return static_cast<int>(n);
    if (n == 0) return static_cast<int>(m);

    // Rolling single-row DP (space O(n))
    std::vector<int> dp(n + 1);
    for (size_t j = 0; j <= n; ++j) dp[j] = static_cast<int>(j);

    for (size_t i = 1; i <= m; ++i) {
        int prev = static_cast<int>(i) - 1;
        dp[0] = static_cast<int>(i);
        for (size_t j = 1; j <= n; ++j) {
            int cur = std::min({
                prev + (a[i-1] == b[j-1] ? 0 : 1),
                dp[j]   + 1,
                dp[j-1] + 1
            });
            prev  = dp[j];
            dp[j] = cur;
        }
    }
    return dp[n];
}

void path_sequence(const VariationGraph& vg,
                   g2g::span<const NodeId> path,
                   std::string& out) {
    for (NodeId n : path) {
        auto sv = vg.seq(n);
        out.append(sv.data(), sv.size());
    }
}

void allele_sequences(const VariationGraph& vg,
                      const Snarl& s,
                      const SnarlTree& tree,
                      std::vector<std::string>& out) {
    // Port of Python Snarl.allele_sequences(vg):
    //   for p in self.allele_paths:
    //       inner = p[1:-1]          # exclude source and sink (shared boundaries)
    //       results.append(path_sequence(vg, inner) if inner else "")
    //
    // The SnarlTree flat store holds full paths (source .. sink inclusive).
    // We skip the first and last node to get inner nodes only.
    const uint32_t n = s.allele_count();
    out.clear();
    out.resize(n);
    for (uint32_t i = 0; i < n; ++i) {
        auto full = tree.allele_path(s.alleles_begin, i);
        out[i].clear();
        // Inner nodes: indices 1 .. size-2  (skip source at 0, sink at end)
        if (full.size() > 2) {
            auto inner = full.subspan(1, full.size() - 2);
            path_sequence(vg, inner, out[i]);
        }
        // If full.size() <= 2, the allele has no inner nodes → empty string
    }
}

} // namespace g2g
