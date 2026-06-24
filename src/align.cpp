/**
 * src/align.cpp
 * =============
 * Graph-to-graph alignment pipeline.
 *
 * P3-6 (allele_sequences wiring) — COMPLETE:
 *   score_snarl_pair() now calls allele_sequences() for both snarls,
 *   runs hungarian_match(), and populates w_struct, allele_matches,
 *   a_only, b_only on the SnarlAlignment.
 *
 * P3-7 (build_deltas) — COMPLETE:
 *   build_deltas() populates AlignmentResult::a_only_alleles and
 *   b_only_alleles from the completed snarl alignments.
 */

#include <g2g/align.hpp>
#include <g2g/anchor.hpp>
#include <g2g/sequence.hpp>
#include <g2g/matching.hpp>

#include <algorithm>
#include <chrono>
#include <deque>
#include <string>
#include <unordered_map>
#include <utility>

namespace g2g {

// ── Snarl correspondence ──────────────────────────────────────────────────────

std::vector<std::pair<uint32_t, uint32_t>> correspond_snarls(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb) {

    // Key a snarl by the bytes of its source and sink sequences.  '\x01' is a
    // separator that cannot occur in DNA, so distinct (source, sink) pairs map
    // to distinct keys.
    auto boundary_key = [](const VariationGraph& vg, const Snarl& s) {
        std::string_view src = vg.seq(s.source);
        std::string_view snk = vg.seq(s.sink);
        std::string key;
        key.reserve(src.size() + snk.size() + 1);
        key.append(src.data(), src.size());
        key.push_back('\x01');
        key.append(snk.data(), snk.size());
        return key;
    };

    // Index B's snarls by boundary key, preserving snarl order so that
    // colliding keys are consumed positionally (FIFO).
    std::unordered_map<std::string, std::deque<uint32_t>> b_by_boundary;
    for (uint32_t ib = 0; ib < tb.snarls.size(); ++ib)
        b_by_boundary[boundary_key(gb, tb.snarls[ib])].push_back(ib);

    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    pairs.reserve(std::min(ta.snarls.size(), tb.snarls.size()));
    for (uint32_t ia = 0; ia < ta.snarls.size(); ++ia) {
        auto it = b_by_boundary.find(boundary_key(ga, ta.snarls[ia]));
        if (it != b_by_boundary.end() && !it->second.empty()) {
            pairs.emplace_back(ia, it->second.front());
            it->second.pop_front();
        }
        // No matching boundary in B → snarl is left unpaired.
    }
    return pairs;
}

AlignmentResult GraphAligner::align(const VariationGraph& ga,
                                     const VariationGraph& gb) const {
    AlignmentResult result;

    using Clock = std::chrono::steady_clock;
    auto secs = [](Clock::time_point a, Clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    };

    // Step 1 — anchor detection
    auto t_idx0 = Clock::now();
    result.anchors = find_anchors(ga, gb, params_.anchor);
    auto t_idx1 = Clock::now();

    // Step 2 — snarl decomposition (allele paths now stored in SnarlTree)
    SnarlTree ta = SnarlDecomposer(ga).decompose();
    SnarlTree tb = SnarlDecomposer(gb).decompose();
    auto t_dec1 = Clock::now();

    // Step 3 — snarl alignment
    result.snarl_alignments = align_snarls(ga, gb, ta, tb, result.anchors);
    auto t_dp1 = Clock::now();

    phase_times_.index_seconds     = secs(t_idx0, t_idx1);
    phase_times_.decompose_seconds = secs(t_idx1, t_dec1);
    phase_times_.dp_seconds        = secs(t_dec1, t_dp1);

    // Step 4 — total score: anchor contributions + snarl scores
    for (const auto& a : result.anchors)
        result.total_score += params_.alpha * a.sim;
    for (const auto& sa : result.snarl_alignments)
        result.total_score += sa.score;

    // Step 5 — delta graph (P3-7)
    build_deltas(ga, gb, result, ta, tb);

    return result;
}

std::vector<SnarlAlignment> GraphAligner::align_snarls(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree& ta,
    const SnarlTree& tb,
    const std::vector<AnchorMatch>&) const {

    std::vector<SnarlAlignment> results;
    std::vector<bool> used_b(tb.snarls.size(), false);

    for (uint32_t ia = 0; ia < ta.snarls.size(); ++ia) {
        float best_score = 0.f;
        int   best_ib    = -1;

        for (uint32_t ib = 0; ib < tb.snarls.size(); ++ib) {
            if (used_b[ib]) continue;
            SnarlAlignment aln = score_snarl_pair(ga, gb, ta, tb,
                                                   ia, ib);
            if (aln.score > best_score) {
                best_score = aln.score;
                best_ib    = static_cast<int>(ib);
            }
        }
        if (best_ib >= 0) {
            results.push_back(score_snarl_pair(ga, gb, ta, tb,
                                                ia, static_cast<uint32_t>(best_ib)));
            used_b[static_cast<size_t>(best_ib)] = true;
        }
    }
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

    // W_topo: 1.0 for simple bubbles (all current snarls qualify)
    aln.w_topo = params_.beta * 1.f;

    // W_struct (P3-6): bipartite allele matching via LAPJV ────────────
    // Port of Python _score_snarl_pair():
    //   seqs_a = sa.allele_sequences(ga)
    //   seqs_b = sb.allele_sequences(gb)
    //   struct_raw, allele_matches = bipartite_allele_match(seqs_a, seqs_b)
    //   w_struct = gamma * struct_raw
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

    aln.w_struct      = params_.gamma * mr.score;
    aln.allele_matches = mr.matches;
    aln.a_only        = mr.a_unmatched;
    aln.b_only        = mr.b_unmatched;

    aln.score = aln.w_seq + aln.w_topo + aln.w_struct;
    return aln;
}

void GraphAligner::build_deltas(
    const VariationGraph& ga,
    const VariationGraph& gb,
    AlignmentResult& result,
    const SnarlTree& ta,
    const SnarlTree& tb) const {

    // Port of Python _build_deltas():
    //   for aln in result.snarl_alignments:
    //       for idx in aln.a_only_alleles:
    //           result.a_only_alleles.append((sa.id, sa.allele_paths[idx]))
    //       for idx in aln.b_only_alleles:
    //           result.b_only_alleles.append((sb.id, sb.allele_paths[idx]))
    //
    // For each unmatched allele, we store:
    //   - snarl_id: "<source_name>><sink_name>"
    //   - path:     the full NodeId path (source..sink inclusive)
    //   - sequence: concatenated inner-node sequences (same as allele_sequences)

    auto node_name = [](const VariationGraph& vg, NodeId n) -> std::string {
        if (n < vg.node_names.size()) return vg.node_names[n];
        return std::to_string(n);
    };

    auto make_snarl_id = [&](const VariationGraph& vg, const Snarl& s) -> std::string {
        return node_name(vg, s.source) + ">" + node_name(vg, s.sink);
    };

    auto inner_seq = [](const VariationGraph& vg,
                         std::span<const NodeId> full_path) -> std::string {
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

    for (const auto& aln : result.snarl_alignments) {
        if (aln.snarl_a_idx >= ta.snarls.size()) continue;
        if (aln.snarl_b_idx >= tb.snarls.size()) continue;

        const Snarl& sa = ta.snarls[aln.snarl_a_idx];
        const Snarl& sb = tb.snarls[aln.snarl_b_idx];

        // A-only alleles: in A, not matched in B
        for (uint16_t idx : aln.a_only) {
            auto path = ta.allele_path(sa.alleles_begin, idx);
            DeltaAllele da;
            da.snarl_id = make_snarl_id(ga, sa);
            da.path     = std::vector<NodeId>(path.begin(), path.end());
            da.sequence = inner_seq(ga, path);
            result.a_only_alleles.push_back(std::move(da));
        }

        // B-only alleles: in B, not matched in A
        for (uint16_t idx : aln.b_only) {
            auto path = tb.allele_path(sb.alleles_begin, idx);
            DeltaAllele db;
            db.snarl_id = make_snarl_id(gb, sb);
            db.path     = std::vector<NodeId>(path.begin(), path.end());
            db.sequence = inner_seq(gb, path);
            result.b_only_alleles.push_back(std::move(db));
        }
    }
}

} // namespace g2g
