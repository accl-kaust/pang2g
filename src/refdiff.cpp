/**
 * src/refdiff.cpp
 * ===============
 * Reference-coordinate graph-to-graph diff (see refdiff.hpp).
 *
 * Each snarl bubble is decomposed into atomic variants (SNPs / indels) by
 * aligning every allele to the snarl's reference allele and keying each
 * difference by its ABSOLUTE reference position.  Because that position is a
 * property of the reference coordinate — not of how a graph chops its nodes or
 * how many adjacent variants a bubble happens to span — a merged bubble in one
 * graph decomposes into the same atomic variants as the separate bubbles in the
 * other, so corresponding variants match.
 */

#include <g2g/refdiff.hpp>
#include <g2g/sequence.hpp>

#include <algorithm>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace g2g {
namespace {

// ── PanSN reference-path resolution ───────────────────────────────────────────

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
    std::size_t colon = contig.find(':');
    p.contig = (colon == std::string::npos) ? contig : contig.substr(0, colon);
    return p;
}

struct RefPath { uint32_t path_idx; std::string contig; };

std::vector<RefPath> resolve_ref_paths(const VariationGraph& g,
                                       const std::string& query) {
    std::vector<RefPath> exact, sample, contig, prefix;
    for (uint32_t p = 0; p < g.num_paths(); ++p) {
        const std::string& name = g.path_names[p];
        PanSN ps = parse_pansn(name);
        RefPath rp{p, ps.contig};
        if (name == query)                                 exact.push_back(rp);
        else if (!ps.sample.empty() && ps.sample == query) sample.push_back(rp);
        else if (ps.contig == query)                       contig.push_back(rp);
        else if (name.rfind(query, 0) == 0)                prefix.push_back(rp);
    }
    if (!exact.empty())  return exact;
    if (!sample.empty()) return sample;
    if (!contig.empty()) return contig;
    return prefix;
}

// ── Allele decomposition into atomic variants ─────────────────────────────────

struct AtomicVar {
    uint32_t    roff;   // 0-based offset within the reference allele
    std::string ref;    // reference bases (empty = insertion)
    std::string alt;    // alternate bases (empty = deletion)
};

// Global (edit-distance) alignment of `alt` against `ref`, walked into atomic
// variants.  Deterministic tie-breaking (diagonal, then deletion, then
// insertion) so the same pair always yields the same decomposition.  Very long
// alleles (large SVs) are left whole to bound cost.
std::vector<AtomicVar> decompose_allele(const std::string& R,
                                        const std::string& A) {
    const std::size_t n = R.size(), m = A.size();
    if (R == A) return {};
    if ((n + 1) * (m + 1) > 9'000'000) {            // huge SV → keep whole
        return {AtomicVar{0, R, A}};
    }

    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
    for (std::size_t i = 0; i <= n; ++i) dp[i][0] = static_cast<int>(i);
    for (std::size_t j = 0; j <= m; ++j) dp[0][j] = static_cast<int>(j);
    for (std::size_t i = 1; i <= n; ++i)
        for (std::size_t j = 1; j <= m; ++j) {
            int sub = dp[i - 1][j - 1] + (R[i - 1] == A[j - 1] ? 0 : 1);
            dp[i][j] = std::min({sub, dp[i - 1][j] + 1, dp[i][j - 1] + 1});
        }

    // Traceback → op string (forward order after reversing).
    std::string ops;
    std::size_t i = n, j = m;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 &&
            dp[i][j] == dp[i - 1][j - 1] + (R[i - 1] == A[j - 1] ? 0 : 1)) {
            ops.push_back('M'); --i; --j;
        } else if (i > 0 && dp[i][j] == dp[i - 1][j] + 1) {
            ops.push_back('D'); --i;            // ref base deleted in alt
        } else {
            ops.push_back('I'); --j;            // base inserted in alt
        }
    }
    std::reverse(ops.begin(), ops.end());

    std::vector<AtomicVar> out;
    std::size_t rp = 0, ap = 0, k = 0;
    while (k < ops.size()) {
        if (ops[k] == 'M') {
            if (R[rp] != A[ap])
                out.push_back({static_cast<uint32_t>(rp),
                               std::string(1, R[rp]), std::string(1, A[ap])});
            ++rp; ++ap; ++k;
        } else if (ops[k] == 'D') {
            std::size_t start = rp; std::string del;
            while (k < ops.size() && ops[k] == 'D') { del.push_back(R[rp]); ++rp; ++k; }
            out.push_back({static_cast<uint32_t>(start), std::move(del), ""});
        } else {
            std::size_t start = rp; std::string ins;
            while (k < ops.size() && ops[k] == 'I') { ins.push_back(A[ap]); ++ap; ++k; }
            out.push_back({static_cast<uint32_t>(start), "", std::move(ins)});
        }
    }
    return out;
}

// (abs reference pos, ref bases, alt bases)
using VarKey = std::tuple<uint64_t, std::string, std::string>;

// Collect every atomic variant along one reference path, mapped to the snarl it
// came from (for reporting matched snarl pairs).
std::map<VarKey, uint32_t> collect_variants(const VariationGraph& g,
                                            const SnarlTree& t,
                                            uint32_t ref_idx) {
    std::map<VarKey, uint32_t> vars;
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

        uint64_t open = off[its->second] + g.seq(s.source).size();

        std::string ref_allele;
        for (uint32_t i = its->second + 1; i < itk->second; ++i) {
            std::string_view sv = g.seq(refnodes[i]);
            ref_allele.append(sv.data(), sv.size());
        }

        std::vector<std::string> alleles;
        allele_sequences(g, s, t, alleles);
        for (const auto& alt : alleles) {
            if (alt == ref_allele) continue;
            for (const AtomicVar& v : decompose_allele(ref_allele, alt))
                vars.emplace(VarKey{open + v.roff, v.ref, v.alt}, si);
        }
    }
    return vars;
}

DeltaAllele make_delta(const std::string& chrom, const VarKey& k) {
    DeltaAllele d;
    d.chrom    = chrom;
    d.pos      = static_cast<int64_t>(std::get<0>(k));
    d.ref      = std::get<1>(k);
    d.sequence = std::get<2>(k);
    d.snarl_id = chrom + ":" + std::to_string(d.pos);
    return d;
}

void diff_contig(const VariationGraph& ga, const SnarlTree& ta, uint32_t a_idx,
                 const VariationGraph& gb, const SnarlTree& tb, uint32_t b_idx,
                 const std::string& chrom, AlignmentResult& res,
                 std::set<std::pair<uint32_t, uint32_t>>& matched_snarls) {
    std::map<VarKey, uint32_t> A = collect_variants(ga, ta, a_idx);
    std::map<VarKey, uint32_t> B = collect_variants(gb, tb, b_idx);

    for (const auto& [key, a_snarl] : A) {
        auto itb = B.find(key);
        if (itb != B.end()) {
            matched_snarls.emplace(a_snarl, itb->second);
            res.total_score += 1.f;          // shared atomic variant
        } else {
            res.a_only_alleles.push_back(make_delta(chrom, key));
        }
    }
    for (const auto& [key, b_snarl] : B) {
        (void)b_snarl;
        if (!A.count(key))
            res.b_only_alleles.push_back(make_delta(chrom, key));
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

    std::unordered_map<std::string, uint32_t> b_by_contig;
    for (const auto& r : rb) b_by_contig.emplace(r.contig, r.path_idx);

    AlignmentResult res;
    std::set<std::pair<uint32_t, uint32_t>> matched_snarls;
    bool any = false;
    for (const auto& r : ra) {
        auto it = b_by_contig.find(r.contig);
        if (it == b_by_contig.end()) continue;
        diff_contig(ga, ta, r.path_idx, gb, tb, it->second, r.contig,
                    res, matched_snarls);
        any = true;
    }
    if (!any) return std::nullopt;

    // Report each matched (A-snarl, B-snarl) pair once, for the GAF.
    for (const auto& [as, bs] : matched_snarls) {
        SnarlAlignment al{};
        al.snarl_a_idx = as;
        al.snarl_b_idx = bs;
        res.snarl_alignments.push_back(al);
    }
    return res;
}

}  // namespace g2g
