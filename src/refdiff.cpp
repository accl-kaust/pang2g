/**
 * src/refdiff.cpp
 * ===============
 * Reference-coordinate graph-to-graph diff (see refdiff.hpp).
 */

#include <g2g/refdiff.hpp>
#include <g2g/sequence.hpp>

#include <algorithm>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace g2g {
namespace {

// One variant site, keyed by the reference position at which the bubble opens.
struct Site {
    uint32_t                 snarl_idx;     // index into the SnarlTree
    std::string              ref_allele;    // inner sequence along the reference
    std::vector<std::string> alts;          // non-reference inner sequences
};

// Collect variant sites of one graph, keyed by reference offset.  Returns false
// if the graph does not embed `ref_path`.
bool collect_sites(const VariationGraph& g,
                   const SnarlTree& t,
                   const std::string& ref_path,
                   std::map<uint64_t, Site>& out) {
    int ref_idx = -1;
    for (uint32_t p = 0; p < g.num_paths(); ++p)
        if (g.path_names[p] == ref_path) { ref_idx = static_cast<int>(p); break; }
    if (ref_idx < 0) return false;

    auto refnodes = g.path_nodes(static_cast<uint32_t>(ref_idx));
    std::unordered_map<NodeId, uint32_t> idx_in_ref;  // node -> position in ref path
    std::vector<uint64_t> off(refnodes.size());       // base offset of each ref node
    uint64_t cum = 0;
    for (uint32_t i = 0; i < refnodes.size(); ++i) {
        idx_in_ref[refnodes[i]] = i;
        off[i] = cum;
        cum += g.seq(refnodes[i]).size();
    }

    for (uint32_t si = 0; si < t.snarls.size(); ++si) {
        const Snarl& s = t.snarls[si];
        auto its = idx_in_ref.find(s.source);
        auto itk = idx_in_ref.find(s.sink);
        if (its == idx_in_ref.end() || itk == idx_in_ref.end())
            continue;  // snarl not anchored on the reference path
        if (itk->second <= its->second)
            continue;  // sink must follow source along the reference

        uint64_t pos = off[its->second] + g.seq(s.source).size();

        Site site;
        site.snarl_idx = si;
        for (uint32_t i = its->second + 1; i < itk->second; ++i) {
            std::string_view sv = g.seq(refnodes[i]);
            site.ref_allele.append(sv.data(), sv.size());
        }

        std::vector<std::string> alleles;
        allele_sequences(g, s, t, alleles);
        for (auto& a : alleles)
            if (a != site.ref_allele)
                site.alts.push_back(std::move(a));

        // Keep the first snarl that opens at a given reference position.
        out.emplace(pos, std::move(site));
    }
    return true;
}

DeltaAllele make_delta(const std::string& chrom, uint64_t pos,
                       const std::string& seq) {
    DeltaAllele d;
    d.chrom    = chrom;
    d.pos      = static_cast<int64_t>(pos);
    d.snarl_id = chrom + ":" + std::to_string(pos);
    d.sequence = seq;
    return d;
}

}  // namespace

std::optional<AlignmentResult> reference_diff(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree&      ta,
    const SnarlTree&      tb,
    const std::string&    ref_path) {

    std::map<uint64_t, Site> A, B;
    if (!collect_sites(ga, ta, ref_path, A)) return std::nullopt;
    if (!collect_sites(gb, tb, ref_path, B)) return std::nullopt;

    AlignmentResult res;

    // Walk the union of reference positions.
    auto ia = A.begin();
    auto ib = B.begin();
    while (ia != A.end() || ib != B.end()) {
        bool has_a = ia != A.end();
        bool has_b = ib != B.end();

        uint64_t pos;
        const Site* sa = nullptr;
        const Site* sb = nullptr;
        if (has_a && (!has_b || ia->first <= ib->first)) {
            pos = ia->first;
            sa  = &ia->second;
            if (has_b && ib->first == pos) { sb = &ib->second; ++ib; }
            ++ia;
        } else {
            pos = ib->first;
            sb  = &ib->second;
            ++ib;
        }

        // Exact-sequence allele sets at this site.
        auto to_set = [](const Site* s) {
            std::vector<std::string> v;
            if (s) v = s->alts;
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
            return v;
        };
        std::vector<std::string> av = to_set(sa);
        std::vector<std::string> bv = to_set(sb);

        std::vector<std::string> shared, a_only, b_only;
        std::set_intersection(av.begin(), av.end(), bv.begin(), bv.end(),
                              std::back_inserter(shared));
        std::set_difference(av.begin(), av.end(), bv.begin(), bv.end(),
                            std::back_inserter(a_only));
        std::set_difference(bv.begin(), bv.end(), av.begin(), av.end(),
                            std::back_inserter(b_only));

        for (const auto& s : a_only)
            res.a_only_alleles.push_back(make_delta(ref_path, pos, s));
        for (const auto& s : b_only)
            res.b_only_alleles.push_back(make_delta(ref_path, pos, s));

        // Record a snarl alignment for sites present (and matched) in both, so
        // the GAF/score reflect corresponding variant sites.
        if (sa && sb) {
            SnarlAlignment al{};
            al.snarl_a_idx = sa->snarl_idx;
            al.snarl_b_idx = sb->snarl_idx;
            al.w_struct    = static_cast<float>(shared.size());
            al.score       = al.w_struct;
            res.snarl_alignments.push_back(std::move(al));
            res.total_score += static_cast<float>(shared.size());
        }
    }

    return res;
}

}  // namespace g2g
