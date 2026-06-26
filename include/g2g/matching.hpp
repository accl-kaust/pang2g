#pragma once
/**
 * include/g2g/matching.hpp
 * ========================
 * Bipartite allele matching — C++ replacement for bipartite_allele_match().
 *
 * Phase 2 finding: brute-force is 1,491× slower at k=7 vs k=2; infeasible
 * at k≥8.  This module provides:
 *
 *   hungarian_match()  — O(k^3) via Jonker-Volgenant LAPJV algorithm.
 *                        Mandatory for k > 6.  This is the Phase 3
 *                        highest-priority task (§3.3, item 4).
 *
 *   greedy_match()     — O(k^2 log k) approximation; (1-1/e) ratio.
 *                        Fallback for k > 30 (Phase 2: rare at HPRC scale
 *                        but present in complex SV regions).
 *
 *   brute_match()      — Exact brute-force; only for k ≤ 6 unit tests
 *                        (correctness oracle against the Python prototype).
 *
 * All three share the same MatchResult type so the aligner can swap them.
 */

#include <g2g/types.hpp>

#include <cstdint>
#include <vector>
#include <g2g/span.hpp>   // portable g2g::span (falls back on GCC 9)

namespace g2g {

struct MatchResult {
    float                    score;          // W_struct (before gamma scaling)
    std::vector<AlleleMatch> matches;        // (i, j, sim) pairs
    std::vector<uint16_t>    a_unmatched;   // allele indices in A with no match
    std::vector<uint16_t>    b_unmatched;   // allele indices in B with no match
};

/**
 * Policy constants (mirror Python defaults).
 */
struct MatchParams {
    float gap_cost     = 0.5f;
    float min_sim      = 0.75f;
    int   lapjv_limit  = 30;   // above this k, switch to greedy
};

/**
 * Primary entry point.
 * Dispatches to hungarian_match or greedy_match based on k and params.
 * seqs_a / seqs_b: one string_view per allele (inner-node sequence).
 */
MatchResult allele_match(g2g::span<const std::string_view> seqs_a,
                         g2g::span<const std::string_view> seqs_b,
                         const MatchParams& params = {});

/**
 * Hungarian (LAPJV) exact O(k^3) matching.
 * Builds the k_a × k_b similarity matrix then solves the assignment problem.
 * Pairs with sim < params.min_sim are treated as if their cost is -∞.
 */
MatchResult hungarian_match(g2g::span<const std::string_view> seqs_a,
                             g2g::span<const std::string_view> seqs_b,
                             const MatchParams& params = {});

/**
 * Greedy (1-1/e) approximation.
 * Sort all (i,j) pairs by sim descending; greedily assign while neither
 * index is used and sim ≥ min_sim.
 */
MatchResult greedy_match(g2g::span<const std::string_view> seqs_a,
                          g2g::span<const std::string_view> seqs_b,
                          const MatchParams& params = {});

/**
 * Brute-force — exact, for k ≤ 6 only.
 * Used as a test oracle to verify hungarian_match output.
 * Mirrors bipartite_allele_match() from align.py exactly.
 */
MatchResult brute_match(g2g::span<const std::string_view> seqs_a,
                         g2g::span<const std::string_view> seqs_b,
                         const MatchParams& params = {});

} // namespace g2g
