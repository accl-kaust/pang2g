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

// PanSN path naming: "<sample>#<hap>#<contig>[:<start>-<end>]".  Plain P-line
// names (e.g. "chr1") have no sample; we treat the whole name as the contig.
struct PanSN { std::string sample; std::string contig; };

PanSN parse_pansn(const std::string& name) {
    PanSN p;
    std::size_t h1 = name.find('#');
    std::size_t h2 = (h1 == std::string::npos) ? std::string::npos
                                               : name.find('#', h1 + 1);
    std::string contig;
    if (h1 != std::string::npos && h2 != std::string::npos) {
        p.sample = name.substr(0, h1);
        contig   = name.substr(h2 + 1);
    } else {
        contig = name;
    }
    std::size_t colon = contig.find(':');   // strip any :start-end range
    p.contig = (colon == std::string::npos) ? contig : contig.substr(0, colon);
    return p;
}

// A reference path resolved from a --ref query, tagged with its display contig.
struct RefPath { uint32_t path_idx; std::string contig; };

// Resolve a --ref query to reference path(s).  Priority: exact name > sample
// (PanSN component, e.g. "GRCh38") > contig (e.g. "chr1") > name prefix.  The
// reference is normally identified by its sample, which is unique per contig;
// a query may resolve to several paths when the reference spans contigs.
std::vector<RefPath> resolve_ref_paths(const VariationGraph& g,
                                       const std::string& query) {
    std::vector<RefPath> exact, sample, contig, prefix;
    for (uint32_t p = 0; p < g.num_paths(); ++p) {
        const std::string& name = g.path_names[p];
        PanSN ps = parse_pansn(name);
        RefPath rp{p, ps.contig};
        if (name == query)                          exact.push_back(rp);
        else if (!ps.sample.empty() && ps.sample == query) sample.push_back(rp);
        else if (ps.contig == query)                contig.push_back(rp);
        else if (name.rfind(query, 0) == 0)         prefix.push_back(rp);
    }
    if (!exact.empty())  return exact;
    if (!sample.empty()) return sample;
    if (!contig.empty()) return contig;
    return prefix;
}

// One variant site, keyed by the reference position at which the bubble opens.
struct Site {
    uint32_t                 snarl_idx;
    std::string              ref_allele;
    std::vector<std::string> alts;
};

// Collect variant sites along one reference path, keyed by reference offset.
std::map<uint64_t, Site> collect_sites(const VariationGraph& g,
                                       const SnarlTree& t,
                                       uint32_t ref_idx) {
    std::map<uint64_t, Site> out;
    auto refnodes = g.path_nodes(ref_idx);
    std::unordered_map<NodeId, uint32_t> idx_in_ref;
    std::vector<uint64_t> off(refnodes.size());
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
        if (its == idx_in_ref.end() || itk == idx_in_ref.end()) continue;
        if (itk->second <= its->second) continue;

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
            if (a != site.ref_allele) site.alts.push_back(std::move(a));

        out.emplace(pos, std::move(site));
    }
    return out;
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

// Diff one shared reference contig (path indices into A and B) into `res`.
void diff_contig(const VariationGraph& ga, const SnarlTree& ta, uint32_t a_idx,
                 const VariationGraph& gb, const SnarlTree& tb, uint32_t b_idx,
                 const std::string& chrom, AlignmentResult& res) {
    std::map<uint64_t, Site> A = collect_sites(ga, ta, a_idx);
    std::map<uint64_t, Site> B = collect_sites(gb, tb, b_idx);

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

        auto to_set = [](const Site* s) {
            std::vector<std::string> v;
            if (s) v = s->alts;
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
            return v;
        };
        std::vector<std::string> av = to_set(sa), bv = to_set(sb);

        std::vector<std::string> shared, a_only, b_only;
        std::set_intersection(av.begin(), av.end(), bv.begin(), bv.end(),
                              std::back_inserter(shared));
        std::set_difference(av.begin(), av.end(), bv.begin(), bv.end(),
                            std::back_inserter(a_only));
        std::set_difference(bv.begin(), bv.end(), av.begin(), av.end(),
                            std::back_inserter(b_only));

        for (const auto& s : a_only)
            res.a_only_alleles.push_back(make_delta(chrom, pos, s));
        for (const auto& s : b_only)
            res.b_only_alleles.push_back(make_delta(chrom, pos, s));

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
}

}  // namespace

std::optional<AlignmentResult> reference_diff(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree&      ta,
    const SnarlTree&      tb,
    const std::string&    ref_path) {

    std::vector<RefPath> ra = resolve_ref_paths(ga, ref_path);
    std::vector<RefPath> rb = resolve_ref_paths(gb, ref_path);
    if (ra.empty() || rb.empty()) return std::nullopt;

    // Match reference contigs between the two graphs by name and diff each.
    std::unordered_map<std::string, uint32_t> b_by_contig;
    for (const auto& r : rb) b_by_contig.emplace(r.contig, r.path_idx);

    AlignmentResult res;
    bool any = false;
    for (const auto& r : ra) {
        auto it = b_by_contig.find(r.contig);
        if (it == b_by_contig.end()) continue;   // contig present in A only
        diff_contig(ga, ta, r.path_idx, gb, tb, it->second, r.contig, res);
        any = true;
    }
    if (!any) return std::nullopt;   // no shared reference contig
    return res;
}

}  // namespace g2g
