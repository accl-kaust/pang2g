/**
 * src/align.cpp
 * =============
 * Graph-to-graph alignment pipeline (sequential driver) and the two snarl
 * pairing rules shared with the parallel runtime.
 *
 * Snarl pairing
 * -------------
 *  * correspond_snarls()          - boundary-key rule (default).  Two snarls
 *    correspond when the byte sequences of their source and sink nodes are
 *    identical; repeated keys are consumed positionally (FIFO).  O(|TA|+|TB|).
 *    Its failure mode: a snarl present in only one graph whose boundary key
 *    repeats shifts every later snarl with that key by one position.
 *  * correspond_snarls_interval() - anchor-interval rule (--pairing interval).
 *    The anchors (backbone nodes shared by both graphs) are chained in
 *    topological order; between two consecutive anchors the snarls of A and B
 *    are paired by an order-preserving dynamic program (a Needleman-Wunsch over
 *    the two snarl sequences) whose match score is the fraction of identical
 *    allele sequences plus the boundary similarity.  Errors cannot cross an
 *    anchor, and inside an interval the pairing is optimal for that score.
 *    Cost is sum_i |A_i|*|B_i| over intervals, near-linear when anchors are
 *    dense; an interval larger than kMaxDpCells falls back to the boundary rule.
 *
 * Delta construction (build_deltas) emits (a) the alleles of paired snarls left
 * unmatched by LAPJV and (b) every allele of a snarl that found no partner, so
 * the delta is exactly the set of alleles of one graph unmatched by the
 * alignment (Definition 2 in the paper).
 */

#include <g2g/align.hpp>
#include <g2g/anchor.hpp>
#include <g2g/decompose.hpp>
#include <g2g/sequence.hpp>
#include <g2g/matching.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <deque>
#include <functional>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace g2g {

// -- Boundary-key pairing (default) --------------------------------------------

namespace {

std::string boundary_key(const VariationGraph& vg, const Snarl& s) {
    // '\x01' cannot occur in DNA, so distinct (source, sink) pairs map to
    // distinct keys.
    std::string_view src = vg.seq(s.source);
    std::string_view snk = vg.seq(s.sink);
    std::string key;
    key.reserve(src.size() + snk.size() + 1);
    key.append(src.data(), src.size());
    key.push_back('\x01');
    key.append(snk.data(), snk.size());
    return key;
}

// Boundary-key FIFO pairing restricted to the given snarl subsets (indices).
void pair_by_boundary(const VariationGraph& ga, const VariationGraph& gb,
                      const SnarlTree& ta, const SnarlTree& tb,
                      const std::vector<uint32_t>& A,
                      const std::vector<uint32_t>& B,
                      std::vector<std::pair<uint32_t, uint32_t>>& pairs) {
    std::unordered_map<std::string, std::deque<uint32_t>> b_by_boundary;
    for (uint32_t ib : B)
        b_by_boundary[boundary_key(gb, tb.snarls[ib])].push_back(ib);
    for (uint32_t ia : A) {
        auto it = b_by_boundary.find(boundary_key(ga, ta.snarls[ia]));
        if (it != b_by_boundary.end() && !it->second.empty()) {
            pairs.emplace_back(ia, it->second.front());
            it->second.pop_front();
        }
    }
}

}  // namespace

std::vector<std::pair<uint32_t, uint32_t>> correspond_snarls(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb) {

    std::vector<uint32_t> A(ta.snarls.size()), B(tb.snarls.size());
    std::iota(A.begin(), A.end(), 0u);
    std::iota(B.begin(), B.end(), 0u);
    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    pairs.reserve(std::min(A.size(), B.size()));
    pair_by_boundary(ga, gb, ta, tb, A, B, pairs);
    return pairs;
}

// -- Anchor-interval, order-preserving pairing ---------------------------------

namespace {

// Largest DP table (cells) attempted for one anchor interval before falling
// back to the boundary-key rule.  16M cells = 64 MiB of float scores.
constexpr uint64_t kMaxDpCells = 1ull << 24;

std::vector<uint32_t> topo_positions(const VariationGraph& g) {
    std::vector<NodeId> order = topological_sort(g);
    std::vector<uint32_t> pos(g.num_nodes(), 0);
    for (uint32_t i = 0; i < order.size(); ++i) pos[order[i]] = i;
    return pos;
}

// Sorted 64-bit hashes of each snarl's allele sequences (inner nodes only).
std::vector<std::vector<uint64_t>> allele_hashes(const VariationGraph& g,
                                                 const SnarlTree& t) {
    std::vector<std::vector<uint64_t>> out(t.snarls.size());
    std::vector<std::string> seqs;
    std::hash<std::string_view> h;
    for (uint32_t i = 0; i < t.snarls.size(); ++i) {
        allele_sequences(g, t.snarls[i], t, seqs);
        auto& v = out[i];
        v.clear();
        v.reserve(seqs.size());
        for (const auto& s : seqs) v.push_back(h(std::string_view(s)));
        std::sort(v.begin(), v.end());
    }
    return out;
}

uint32_t sorted_intersection_size(const std::vector<uint64_t>& a,
                                  const std::vector<uint64_t>& b) {
    uint32_t n = 0;
    std::size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] < b[j]) ++i;
        else if (b[j] < a[i]) ++j;
        else { ++n; ++i; ++j; }
    }
    return n;
}

// Longest chain of anchors increasing in both graphs: sort by position in A,
// then take the longest strictly increasing subsequence of positions in B.
std::vector<std::pair<uint32_t, uint32_t>> anchor_chain(
        const std::vector<AnchorMatch>& anchors,
        const std::vector<uint32_t>& pa, const std::vector<uint32_t>& pb) {
    std::vector<std::pair<uint32_t, uint32_t>> v;
    v.reserve(anchors.size());
    for (const auto& an : anchors) v.emplace_back(pa[an.node_a], pb[an.node_b]);
    std::sort(v.begin(), v.end());
    // Patience LIS with predecessor links (strict in .second).
    std::vector<uint32_t> tails;          // index into v of the tail of each length
    std::vector<int64_t> prev(v.size(), -1);
    for (uint32_t i = 0; i < v.size(); ++i) {
        auto it = std::lower_bound(tails.begin(), tails.end(), v[i].second,
                                   [&](uint32_t idx, uint32_t val) {
                                       return v[idx].second < val;
                                   });
        if (it != tails.begin()) prev[i] = *(it - 1);
        if (it == tails.end()) tails.push_back(i);
        else *it = i;
    }
    std::vector<std::pair<uint32_t, uint32_t>> chain;
    if (tails.empty()) return chain;
    for (int64_t k = tails.back(); k >= 0; k = prev[static_cast<std::size_t>(k)])
        chain.push_back(v[static_cast<std::size_t>(k)]);
    std::reverse(chain.begin(), chain.end());
    return chain;
}

// Order-preserving pairing of the A-snarls `A` with the B-snarls `B` (both in
// topological order) maximising the summed pair score; pairs with score 0 are
// not allowed.  Gaps (unpaired snarls) cost nothing.
template <class ScoreFn>
void pair_interval_dp(const std::vector<uint32_t>& A,
                      const std::vector<uint32_t>& B,
                      ScoreFn&& score,
                      std::vector<std::pair<uint32_t, uint32_t>>& pairs) {
    const std::size_t n = A.size(), m = B.size();
    if (n == 0 || m == 0) return;
    const std::size_t W = m + 1;
    std::vector<float>   dp((n + 1) * W, 0.f);
    std::vector<uint8_t> tr((n + 1) * W, 0);   // 0 = up, 1 = left, 2 = diag
    for (std::size_t i = 1; i <= n; ++i) {
        for (std::size_t j = 1; j <= m; ++j) {
            float best = dp[(i - 1) * W + j]; uint8_t t = 0;
            if (dp[i * W + (j - 1)] > best) { best = dp[i * W + (j - 1)]; t = 1; }
            float s = score(A[i - 1], B[j - 1]);
            if (s > 0.f) {
                float d = dp[(i - 1) * W + (j - 1)] + s;
                if (d >= best) { best = d; t = 2; }
            }
            dp[i * W + j] = best; tr[i * W + j] = t;
        }
    }
    std::size_t i = n, j = m;
    std::vector<std::pair<uint32_t, uint32_t>> local;
    while (i > 0 && j > 0) {
        uint8_t t = tr[i * W + j];
        if (t == 2)      { local.emplace_back(A[i - 1], B[j - 1]); --i; --j; }
        else if (t == 1) { --j; }
        else             { --i; }
    }
    pairs.insert(pairs.end(), local.rbegin(), local.rend());
}

}  // namespace

std::vector<std::pair<uint32_t, uint32_t>> correspond_snarls_interval(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb,
    const std::vector<AnchorMatch>& anchors) {

    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    if (ta.snarls.empty() || tb.snarls.empty()) return pairs;

    const std::vector<uint32_t> pa = topo_positions(ga);
    const std::vector<uint32_t> pb = topo_positions(gb);
    const auto ha = allele_hashes(ga, ta);
    const auto hb = allele_hashes(gb, tb);

    // Snarls in topological order of their source node.
    std::vector<uint32_t> sa(ta.snarls.size()), sb(tb.snarls.size());
    std::iota(sa.begin(), sa.end(), 0u);
    std::iota(sb.begin(), sb.end(), 0u);
    std::sort(sa.begin(), sa.end(), [&](uint32_t x, uint32_t y) {
        return pa[ta.snarls[x].source] < pa[ta.snarls[y].source]; });
    std::sort(sb.begin(), sb.end(), [&](uint32_t x, uint32_t y) {
        return pb[tb.snarls[x].source] < pb[tb.snarls[y].source]; });

    const auto chain = anchor_chain(anchors, pa, pb);

    auto score = [&](uint32_t ia, uint32_t ib) -> float {
        uint32_t shared = sorted_intersection_size(ha[ia], hb[ib]);
        if (shared == 0) return 0.f;        // nothing to match → not a candidate
        float frac = static_cast<float>(shared) /
                     static_cast<float>(std::max(ha[ia].size(), hb[ib].size()));
        const Snarl& s = ta.snarls[ia];
        const Snarl& t = tb.snarls[ib];
        float bsim = 0.5f * (seq_similarity(ga.seq(s.source), gb.seq(t.source)) +
                             seq_similarity(ga.seq(s.sink),   gb.seq(t.sink)));
        return frac + 0.5f * bsim;
    };

    std::size_t ia = 0, ib = 0;
    std::vector<uint32_t> A, B;
    for (std::size_t c = 0; c <= chain.size(); ++c) {
        const uint32_t limA = c < chain.size() ? chain[c].first  : UINT32_MAX;
        const uint32_t limB = c < chain.size() ? chain[c].second : UINT32_MAX;
        A.clear(); B.clear();
        while (ia < sa.size() && pa[ta.snarls[sa[ia]].source] < limA) A.push_back(sa[ia++]);
        while (ib < sb.size() && pb[tb.snarls[sb[ib]].source] < limB) B.push_back(sb[ib++]);
        if (A.empty() || B.empty()) continue;
        if (static_cast<uint64_t>(A.size()) * static_cast<uint64_t>(B.size()) > kMaxDpCells)
            pair_by_boundary(ga, gb, ta, tb, A, B, pairs);   // too few anchors here
        else
            pair_interval_dp(A, B, score, pairs);
    }
    std::sort(pairs.begin(), pairs.end());
    return pairs;
}

std::vector<std::pair<uint32_t, uint32_t>> pair_snarls(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb,
    const std::vector<AnchorMatch>& anchors,
    PairingMode mode) {
    if (mode == PairingMode::Interval)
        return correspond_snarls_interval(ga, gb, ta, tb, anchors);
    return correspond_snarls(ga, gb, ta, tb);
}

// -- Delta construction (shared by both aligners) ------------------------------

void build_deltas(const VariationGraph& ga,
                  const VariationGraph& gb,
                  const SnarlTree& ta,
                  const SnarlTree& tb,
                  AlignmentResult& result) {

    auto node_name = [](const VariationGraph& vg, NodeId n) -> std::string {
        if (n < vg.node_names.size()) return vg.node_names[n];
        return std::to_string(n);
    };
    auto make_snarl_id = [&](const VariationGraph& vg, const Snarl& s) -> std::string {
        return node_name(vg, s.source) + ">" + node_name(vg, s.sink);
    };
    auto inner_seq = [](const VariationGraph& vg,
                        g2g::span<const NodeId> full_path) -> std::string {
        std::string seq;
        if (full_path.size() > 2) {
            auto inner = full_path.subspan(1, full_path.size() - 2);
            for (NodeId n : inner) {
                auto sv = vg.seq(n);
                seq.append(sv.data(), sv.size());
            }
        }
        return seq;
    };
    auto emit = [&](const VariationGraph& vg, const SnarlTree& t, const Snarl& s,
                    uint32_t allele_idx, std::vector<DeltaAllele>& out) {
        auto path = t.allele_path(s.alleles_begin, allele_idx);
        DeltaAllele d;
        d.snarl_id = make_snarl_id(vg, s);
        d.path     = std::vector<NodeId>(path.begin(), path.end());
        d.sequence = inner_seq(vg, path);
        out.push_back(std::move(d));
    };

    std::vector<bool> paired_a(ta.snarls.size(), false);
    std::vector<bool> paired_b(tb.snarls.size(), false);

    // (a) alleles of paired snarls left unmatched by the per-snarl assignment
    for (const auto& aln : result.snarl_alignments) {
        if (aln.snarl_a_idx >= ta.snarls.size()) continue;
        if (aln.snarl_b_idx >= tb.snarls.size()) continue;
        paired_a[aln.snarl_a_idx] = true;
        paired_b[aln.snarl_b_idx] = true;
        const Snarl& sa = ta.snarls[aln.snarl_a_idx];
        const Snarl& sb = tb.snarls[aln.snarl_b_idx];
        for (uint16_t idx : aln.a_only) emit(ga, ta, sa, idx, result.a_only_alleles);
        for (uint16_t idx : aln.b_only) emit(gb, tb, sb, idx, result.b_only_alleles);
    }

    // (b) every allele of a snarl that found no partner
    result.unpaired_snarls_a = 0;
    result.unpaired_snarls_b = 0;
    for (uint32_t i = 0; i < ta.snarls.size(); ++i) {
        if (paired_a[i]) continue;
        ++result.unpaired_snarls_a;
        const Snarl& s = ta.snarls[i];
        for (uint32_t k = 0; k < s.allele_count(); ++k) emit(ga, ta, s, k, result.a_only_alleles);
    }
    for (uint32_t i = 0; i < tb.snarls.size(); ++i) {
        if (paired_b[i]) continue;
        ++result.unpaired_snarls_b;
        const Snarl& s = tb.snarls[i];
        for (uint32_t k = 0; k < s.allele_count(); ++k) emit(gb, tb, s, k, result.b_only_alleles);
    }
}

// -- Sequential aligner --------------------------------------------------------

AlignmentResult GraphAligner::align(const VariationGraph& ga,
                                     const VariationGraph& gb) const {
    AlignmentResult result;

    using Clock = std::chrono::steady_clock;
    auto secs = [](Clock::time_point a, Clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    };

    // Step 1 - anchor detection
    auto t_idx0 = Clock::now();
    result.anchors = find_anchors(ga, gb, params_.anchor);
    auto t_idx1 = Clock::now();

    // Step 2 - snarl decomposition (allele paths stored in the SnarlTree)
    SnarlTree ta = SnarlDecomposer(ga).decompose();
    SnarlTree tb = SnarlDecomposer(gb).decompose();
    auto t_dec1 = Clock::now();

    // Step 3 - snarl pairing + per-pair matching
    result.snarl_alignments = align_snarls(ga, gb, ta, tb, result.anchors);
    auto t_dp1 = Clock::now();

    phase_times_.index_seconds     = secs(t_idx0, t_idx1);
    phase_times_.decompose_seconds = secs(t_idx1, t_dec1);
    phase_times_.dp_seconds        = secs(t_dec1, t_dp1);

    // Step 4 - total score: anchor contributions + snarl scores
    for (const auto& a : result.anchors)
        result.total_score += params_.alpha * a.sim;
    for (const auto& sa : result.snarl_alignments)
        result.total_score += sa.score;

    // Step 5 - delta graphs
    build_deltas(ga, gb, ta, tb, result);

    return result;
}

std::vector<SnarlAlignment> GraphAligner::align_snarls(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb,
    const std::vector<AnchorMatch>& anchors) const {

    // The SAME pairing the parallel runtime uses (ParallelAligner::
    // build_task_dag → pair_snarls), so the sequential and parallel aligners
    // produce identical snarl alignments and total scores at every scale.
    std::vector<SnarlAlignment> results;
    auto corr = pair_snarls(ga, gb, ta, tb, anchors, params_.pairing);
    results.reserve(corr.size());
    for (const auto& [ia, ib] : corr)
        results.push_back(score_snarl_pair(ga, gb, ta, tb, ia, ib));
    return results;
}

SnarlAlignment GraphAligner::score_snarl_pair(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb,
    uint32_t ia,
    uint32_t ib) const {

    const Snarl& sa = ta.snarls[ia];
    const Snarl& sb = tb.snarls[ib];

    SnarlAlignment aln{};
    aln.snarl_a_idx = ia;
    aln.snarl_b_idx = ib;

    // W_seq: average similarity of source and sink boundary nodes
    float src_sim  = seq_similarity(ga.seq(sa.source), gb.seq(sb.source));
    float sink_sim = seq_similarity(ga.seq(sa.sink),   gb.seq(sb.sink));
    aln.w_seq = params_.alpha * (src_sim + sink_sim) / 2.f;

    // W_topo: 1.0 for simple bubbles (all depth-1 snarls qualify)
    aln.w_topo = params_.beta * 1.f;

    // W_struct: bipartite allele matching via LAPJV
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

} // namespace g2g
