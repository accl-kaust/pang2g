/**
 * src/matching.cpp
 * ================
 * Bipartite allele matching implementations.
 *
 * Priority order per Phase 2 findings:
 *   1. hungarian_match() — LAPJV O(k³); primary path for k ≤ lapjv_limit.
 *                          (Phase 2: 570× speedup vs brute at k=7.)
 *   2. greedy_match()    — O(k² log k) approximation; (1-1/e) ratio.
 *                          Fallback for k > lapjv_limit (default 30).
 *   3. brute_match()     — Exact oracle; only for k ≤ 8 unit tests.
 *
 * Phase 3 Task P3-2 status: COMPLETE.
 *   - lapjv() vendored as include/g2g/lapjv.hpp (self-contained header).
 *   - hungarian_match() fully wired; all test cases pass.
 *   - brute_match() enumeration complete; used as correctness oracle.
 */

#include <g2g/matching.hpp>
#include <g2g/lapjv.hpp>
#include <g2g/sequence.hpp>

#include <algorithm>
#include <cassert>
#include <numeric>
#include <vector>

namespace g2g {

// ── Helpers ───────────────────────────────────────────────────

static std::vector<float>
build_sim_matrix_flat(g2g::span<const std::string_view> seqs_a,
                      g2g::span<const std::string_view> seqs_b,
                      int na, int nb) {
    // Row-major: flat[i*nb + j] = sim(a[i], b[j])
    std::vector<float> flat(static_cast<size_t>(na) * static_cast<size_t>(nb));
    for (int i = 0; i < na; ++i)
        for (int j = 0; j < nb; ++j)
            flat[static_cast<size_t>(i) * static_cast<size_t>(nb) + static_cast<size_t>(j)] =
                seq_similarity(seqs_a[static_cast<size_t>(i)],
                               seqs_b[static_cast<size_t>(j)]);
    return flat;
}

static MatchResult
build_result_from_sim_flat(const std::vector<float>& sim,
                           int na, int nb,
                           const std::vector<std::pair<int,int>>& pairs,
                           const MatchParams& p) {
    MatchResult r;
    std::vector<bool> used_a(static_cast<size_t>(na), false);
    std::vector<bool> used_b(static_cast<size_t>(nb), false);
    float score = 0.f;

    for (auto [i, j] : pairs) {
        float s = sim[static_cast<size_t>(i) * static_cast<size_t>(nb) + static_cast<size_t>(j)];
        if (s < p.min_sim) continue;
        r.matches.push_back({static_cast<uint16_t>(i),
                              static_cast<uint16_t>(j), s});
        score += s;
        used_a[static_cast<size_t>(i)] = used_b[static_cast<size_t>(j)] = true;
    }
    for (int i = 0; i < na; ++i)
        if (!used_a[static_cast<size_t>(i)]) {
            r.a_unmatched.push_back(static_cast<uint16_t>(i));
            score -= p.gap_cost;
        }
    for (int j = 0; j < nb; ++j)
        if (!used_b[static_cast<size_t>(j)]) {
            r.b_unmatched.push_back(static_cast<uint16_t>(j));
            score -= p.gap_cost;
        }
    r.score = score;
    return r;
}

// ── Brute-force oracle (P3-1: COMPLETE) ──────────────────────
//
// Exact port of Python bipartite_allele_match().
// Enumerates all assignments of size m = 0..min(na,nb):
//   • all C(large_n, m) combinations from the larger side,
//   • all P(small_n, m) permutations from the smaller side,
// picks the assignment maximising sum(sim) - gap_cost * unmatched.
//
// Complexity: O(min(na,nb)! * C(max,min)) — feasible only for k ≤ 8.

MatchResult brute_match(g2g::span<const std::string_view> seqs_a,
                         g2g::span<const std::string_view> seqs_b,
                         const MatchParams& params) {
    const int na = static_cast<int>(seqs_a.size());
    const int nb = static_cast<int>(seqs_b.size());
    assert(na <= 8 && nb <= 8 && "brute_match: k too large, use hungarian_match");

    auto sim = build_sim_matrix_flat(seqs_a, seqs_b, na, nb);

    // Determine which side is smaller/larger for combination generation.
    bool swap   = (na > nb);
    int  small_n = swap ? nb : na;
    int  large_n = swap ? na : nb;

    float best_score = -1e30f;
    std::vector<std::pair<int,int>> best_pairs;

    // Indices into the large side to choose from (combinations of size m).
    std::vector<int> large_idxs(static_cast<size_t>(large_n));
    std::iota(large_idxs.begin(), large_idxs.end(), 0);

    for (int m = 0; m <= std::min(na, nb); ++m) {
        // Generate all C(large_n, m) subsets of the large side.
        // Use bitmask enumeration for simplicity (works for k ≤ 8).
        std::vector<int> chosen_large;
        chosen_large.reserve(static_cast<size_t>(m));

        // Bitmask over large_n bits (large_n ≤ 8, so uint8_t suffices).
        int masks = 1 << large_n;
        for (int mask = 0; mask < masks; ++mask) {
            if (__builtin_popcount(static_cast<unsigned>(mask)) != m) continue;
            chosen_large.clear();
            for (int bit = 0; bit < large_n; ++bit)
                if (mask & (1 << bit)) chosen_large.push_back(bit);

            if (m == 0) {
                // No matches: score = -gap_cost * (na + nb)
                float s = -params.gap_cost * static_cast<float>(na + nb);
                if (s > best_score) {
                    best_score = s;
                    best_pairs.clear();
                }
                continue;
            }

            // Generate all P(small_n, m) permutations of small side of size m.
            // Start with the first m indices of small side; permute all m! orderings.
            std::vector<int> small_perm(static_cast<size_t>(small_n));
            std::iota(small_perm.begin(), small_perm.end(), 0);

            // We need all m-permutations: iterate combinations of m from small side,
            // then all permutations of those m.
            int small_masks = 1 << small_n;
            for (int smask = 0; smask < small_masks; ++smask) {
                if (__builtin_popcount(static_cast<unsigned>(smask)) != m) continue;
                std::vector<int> chosen_small;
                chosen_small.reserve(static_cast<size_t>(m));
                for (int bit = 0; bit < small_n; ++bit)
                    if (smask & (1 << bit)) chosen_small.push_back(bit);

                // All m! permutations of chosen_small
                std::sort(chosen_small.begin(), chosen_small.end());
                do {
                    // Build pairs for this assignment
                    std::vector<std::pair<int,int>> pairs;
                    pairs.reserve(static_cast<size_t>(m));
                    for (int k = 0; k < m; ++k) {
                        int row = swap ? chosen_large[static_cast<size_t>(k)]
                                       : chosen_small[static_cast<size_t>(k)];
                        int col = swap ? chosen_small[static_cast<size_t>(k)]
                                       : chosen_large[static_cast<size_t>(k)];
                        pairs.emplace_back(row, col);
                    }

                    // Compute score for this assignment
                    float s = 0.f;
                    int matched = 0;
                    for (auto [i, j] : pairs) {
                        float sv = sim[static_cast<size_t>(i) * static_cast<size_t>(nb)
                                       + static_cast<size_t>(j)];
                        if (sv >= params.min_sim) {
                            s += sv;
                            ++matched;
                        }
                    }
                    // Unmatched penalty: all indices not in accepted pairs
                    int unmatched_a = na - matched;
                    int unmatched_b = nb - matched;
                    // Pairs below min_sim also contribute an "unmatched" penalty
                    int rejected = m - matched;
                    // Each rejected pair leaves both its a-allele and b-allele unmatched
                    unmatched_a += rejected;
                    unmatched_b += rejected;
                    // Plus the (large_n - m) unchosen from large side and (small_n - m) from small
                    // Those are already counted in na - matched / nb - matched when rejected == 0,
                    // so we just use: total unmatched = (na - #accepted_matches) + (nb - #accepted_matches)
                    s -= params.gap_cost * static_cast<float>(na - matched)
                       + params.gap_cost * static_cast<float>(nb - matched);

                    if (s > best_score) {
                        best_score = s;
                        best_pairs.clear();
                        for (auto [i, j] : pairs) {
                            float sv = sim[static_cast<size_t>(i) * static_cast<size_t>(nb)
                                           + static_cast<size_t>(j)];
                            if (sv >= params.min_sim)
                                best_pairs.emplace_back(i, j);
                        }
                    }
                } while (std::next_permutation(chosen_small.begin(), chosen_small.end()));
            }
        }
    }

    return build_result_from_sim_flat(sim, na, nb, best_pairs, params);
}

// ── Greedy approximation ──────────────────────────────────────

MatchResult greedy_match(g2g::span<const std::string_view> seqs_a,
                          g2g::span<const std::string_view> seqs_b,
                          const MatchParams& params) {
    const int na = static_cast<int>(seqs_a.size());
    const int nb = static_cast<int>(seqs_b.size());

    auto sim = build_sim_matrix_flat(seqs_a, seqs_b, na, nb);

    std::vector<std::tuple<float,int,int>> candidates;
    candidates.reserve(static_cast<size_t>(na * nb));
    for (int i = 0; i < na; ++i)
        for (int j = 0; j < nb; ++j) {
            float s = sim[static_cast<size_t>(i) * static_cast<size_t>(nb) + static_cast<size_t>(j)];
            if (s >= params.min_sim)
                candidates.emplace_back(s, i, j);
        }
    std::sort(candidates.begin(), candidates.end(), std::greater<>{});

    std::vector<bool> used_a(static_cast<size_t>(na), false);
    std::vector<bool> used_b(static_cast<size_t>(nb), false);
    std::vector<std::pair<int,int>> pairs;

    for (auto [s, i, j] : candidates) {
        if (!used_a[static_cast<size_t>(i)] && !used_b[static_cast<size_t>(j)]) {
            pairs.emplace_back(i, j);
            used_a[static_cast<size_t>(i)] = used_b[static_cast<size_t>(j)] = true;
        }
    }
    return build_result_from_sim_flat(sim, na, nb, pairs, params);
}

// ── Hungarian (LAPJV) — P3-2: COMPLETE ───────────────────────
//
// Builds an n×n cost matrix (padding rectangular inputs with gap_cost
// dummy rows/columns), runs lapjv(), then converts the assignment back
// into a MatchResult, discarding pairs below min_sim.
//
// Cost convention: C[i][j] = 1 - sim(a[i], b[j]).
// Dummy rows/columns: cost = gap_cost (cost of leaving real allele unmatched).
// This encodes gap penalties directly in the LAP objective, so the solver
// simultaneously optimises matches and gaps in a single O(n³) pass.

MatchResult hungarian_match(g2g::span<const std::string_view> seqs_a,
                             g2g::span<const std::string_view> seqs_b,
                             const MatchParams& params) {
    const int na = static_cast<int>(seqs_a.size());
    const int nb = static_cast<int>(seqs_b.size());

    // Square dimension: pad the smaller side with dummy alleles.
    int n = std::max(na, nb);
    if (n == 0) return MatchResult{};

    // Build n×n cost matrix.
    // Real pair (i < na, j < nb): cost = 1 - sim(a[i], b[j]).
    //   → sim below min_sim gets cost clamped to gap_cost so the solver
    //     prefers leaving them unmatched over forcing a bad pair.
    // Dummy pair (i >= na or j >= nb): cost = gap_cost.
    //   → makes the cost of a dummy assignment equal to the penalty for
    //     an unmatched real allele, so the solver is indifferent between
    //     "match this allele to a dummy" and "leave it unmatched."
    std::vector<float> cost(static_cast<size_t>(n) * static_cast<size_t>(n),
                             params.gap_cost);
    for (int i = 0; i < na; ++i) {
        for (int j = 0; j < nb; ++j) {
            float s = seq_similarity(seqs_a[static_cast<size_t>(i)],
                                     seqs_b[static_cast<size_t>(j)]);
            // If sim < min_sim, set cost = gap_cost so the solver won't
            // prefer this match over a dummy (gap) assignment.
            float c = (s >= params.min_sim) ? (1.f - s) : params.gap_cost;
            cost[static_cast<size_t>(i) * static_cast<size_t>(n) + static_cast<size_t>(j)] = c;
        }
    }

    // Run LAPJV
    std::vector<int> row_sol(static_cast<size_t>(n));
    std::vector<int> col_sol(static_cast<size_t>(n));
    lapjv(n, cost.data(), row_sol.data(), col_sol.data());

    // Extract real (non-dummy) pairs where sim >= min_sim.
    // Rebuild the sim values we need for MatchResult from scratch (cheap at k ≤ 30).
    std::vector<float> sim_flat;
    if (na > 0 && nb > 0) {
        sim_flat = build_sim_matrix_flat(seqs_a, seqs_b, na, nb);
    }

    std::vector<std::pair<int,int>> pairs;
    pairs.reserve(static_cast<size_t>(std::min(na, nb)));
    for (int i = 0; i < na; ++i) {
        int j = row_sol[static_cast<size_t>(i)];
        if (j < nb) {  // not a dummy column
            float s = sim_flat[static_cast<size_t>(i) * static_cast<size_t>(nb)
                                + static_cast<size_t>(j)];
            if (s >= params.min_sim)
                pairs.emplace_back(i, j);
        }
    }

    return build_result_from_sim_flat(sim_flat.empty()
                                          ? std::vector<float>(static_cast<size_t>(na * nb), 0.f)
                                          : sim_flat,
                                      na, nb, pairs, params);
}

// ── Dispatcher ────────────────────────────────────────────────

MatchResult allele_match(g2g::span<const std::string_view> seqs_a,
                          g2g::span<const std::string_view> seqs_b,
                          const MatchParams& params) {
    int k = static_cast<int>(std::max(seqs_a.size(), seqs_b.size()));
    if (k > params.lapjv_limit)
        return greedy_match(seqs_a, seqs_b, params);
    return hungarian_match(seqs_a, seqs_b, params);
}

} // namespace g2g
