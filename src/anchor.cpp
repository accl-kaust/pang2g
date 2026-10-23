/**
 * src/anchor.cpp
 * ==============
 * Anchor node detection between two variation graphs.
 *
 * Phase 3 Task P3-5b: MinimizerIndex - COMPLETE.
 *
 * Design (from §2.5.4 of the paper):
 *   Build a minimizer sketch index (sliding-window k-mer sketches, as in
 *   minimap2) over graph A's backbone node sequences.  Anchor candidates
 *   are node pairs sharing ≥1 minimizer.  Refine with full seq_similarity;
 *   keep pairs above min_sim.
 *
 * Hash function: FNV-1a 32-bit over the k-mer characters.
 *   Chosen over std::hash for determinism and good avalanche on DNA alphabets
 *   (only 4 characters → std::hash shows clustering on short strings).
 *
 * Minimizer selection: sliding window of size w over sequence positions.
 *   For each window, the minimizer is the hash-minimum k-mer in the window.
 *   Same convention as minimap2 (Li 2018), which ensures that two identical
 *   sequences always share minimizers regardless of window alignment.
 *
 * Threshold: find_anchors() uses MinimizerIndex when
 *   |backbone_a| * |backbone_b| > kNaiveThreshold (default 10^6).
 *   Below the threshold the naïve O(|VA||VB|) scan is faster due to
 *   lower constant factors.
 *
 * Sensitivity analysis (paper §2.7.3, ablation Q3):
 *   The test suite includes a SensitivityAnalysis test that sweeps
 *   kmer_len ∈ {11,15,21} and window ∈ {5,10,15} on synthetic graphs
 *   with planted anchors and reports recall at each setting.
 *   These numbers feed directly into the ablation figure.
 */

#include <g2g/anchor.hpp>
#include <g2g/sequence.hpp>

#include <algorithm>
#include <cassert>
#include <limits>

namespace g2g {

// -- Constants -------------------------------------------------

static constexpr uint64_t kNaiveThreshold = 1'000'000ULL;

// -- FNV-1a 32-bit hash ----------------------------------------

uint32_t MinimizerIndex::hash_kmer(std::string_view kmer) {
    // FNV-1a: fast, deterministic, good avalanche on short DNA strings.
    uint32_t hash = 2166136261u;
    for (char c : kmer) {
        hash ^= static_cast<uint8_t>(c);
        hash *= 16777619u;
    }
    return hash;
}

// -- Minimizer extraction --------------------------------------

std::vector<uint32_t> MinimizerIndex::minimizers(std::string_view seq) const {
    const int k = params_.kmer_len;
    const int w = params_.window;
    const int n = static_cast<int>(seq.size());

    if (n < k) return {};

    // Compute hash of every k-mer
    int num_kmers = n - k + 1;
    std::vector<uint32_t> hashes(static_cast<size_t>(num_kmers));
    for (int i = 0; i < num_kmers; ++i)
        hashes[static_cast<size_t>(i)] = hash_kmer(seq.substr(
            static_cast<size_t>(i), static_cast<size_t>(k)));

    // Sliding window minimum: for each window [i, i+w), emit the minimum hash.
    // Deduplicate consecutive identical minimizers.
    std::vector<uint32_t> result;
    uint32_t prev = std::numeric_limits<uint32_t>::max();

    int num_windows = num_kmers - w + 1;
    if (num_windows <= 0) {
        // Sequence shorter than one full window: take the global minimum
        uint32_t m = *std::min_element(hashes.begin(), hashes.end());
        return {m};
    }

    for (int i = 0; i < num_windows; ++i) {
        uint32_t m = *std::min_element(
            hashes.begin() + i,
            hashes.begin() + i + w);
        if (m != prev) {
            result.push_back(m);
            prev = m;
        }
    }
    return result;
}

// -- Index construction ----------------------------------------

MinimizerIndex::MinimizerIndex(const VariationGraph& vg, const AnchorParams& p)
    : vg_(vg), params_(p) { build(); }

void MinimizerIndex::build() {
    const uint32_t n = vg_.num_nodes();
    backbone_.clear();
    index_.clear();

    for (NodeId nd = 0; nd < n; ++nd) {
        if (vg_.in_degree(nd) <= 1 && vg_.out_degree(nd) <= 1)
            backbone_.push_back(nd);
    }

    for (NodeId nd : backbone_) {
        auto mins = minimizers(vg_.seq(nd));
        for (uint32_t h : mins)
            index_[h].push_back(nd);
    }

    // Deduplicate per-bucket (a node can produce the same minimizer multiple
    // times from different windows; we want it once per bucket)
    for (auto& [h, ids] : index_) {
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    }
}

// -- Query -----------------------------------------------------

std::vector<NodeId> MinimizerIndex::query(std::string_view query_seq) const {
    auto mins = minimizers(query_seq);

    std::vector<NodeId> candidates;
    for (uint32_t h : mins) {
        auto it = index_.find(h);
        if (it != index_.end())
            candidates.insert(candidates.end(),
                               it->second.begin(), it->second.end());
    }

    // Deduplicate candidates
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()),
                     candidates.end());
    return candidates;
}

// -- find_anchors ----------------------------------------------

std::vector<AnchorMatch> find_anchors(const VariationGraph& ga,
                                       const VariationGraph& gb,
                                       const AnchorParams& params) {
    // Collect backbone nodes
    std::vector<NodeId> backbone_a, backbone_b;
    for (NodeId n = 0; n < ga.num_nodes(); ++n)
        if (ga.in_degree(n) <= 1 && ga.out_degree(n) <= 1)
            backbone_a.push_back(n);
    for (NodeId n = 0; n < gb.num_nodes(); ++n)
        if (gb.in_degree(n) <= 1 && gb.out_degree(n) <= 1)
            backbone_b.push_back(n);

    uint64_t naive_work = static_cast<uint64_t>(backbone_a.size())
                        * static_cast<uint64_t>(backbone_b.size());

    std::vector<AnchorMatch> anchors;
    std::vector<bool> used_b(gb.num_nodes(), false);

    if (naive_work <= kNaiveThreshold
#ifdef G2G_DISABLE_MINIMIZER_INDEX
        || true   // ablation A4: always use naive path
#endif
    ) {
        // -- Naïve O(|backbone_a| * |backbone_b|) scan ---------
        for (NodeId na : backbone_a) {
            float  best_sim = params.min_sim - 1e-9f;
            NodeId best_nb  = kInvalidNode;
            for (NodeId nb : backbone_b) {
                if (used_b[nb]) continue;
                float s = seq_similarity(ga.seq(na), gb.seq(nb));
                if (s > best_sim) { best_sim = s; best_nb = nb; }
            }
            if (best_nb != kInvalidNode) {
                anchors.push_back({na, best_nb, best_sim});
                used_b[best_nb] = true;
            }
        }
    } else {
        // -- MinimizerIndex path --------------------------------
        // Index graph B; query each node in graph A's backbone.
        // Candidate pairs share ≥1 minimizer; refine with full seq_similarity.
        MinimizerIndex idx_b(gb, params);

        // Collect all candidates with their exact similarity
        struct Candidate { NodeId na, nb; float sim; };
        std::vector<Candidate> cands;
        cands.reserve(backbone_a.size() * 4);  // rough estimate

        for (NodeId na : backbone_a) {
            auto candidates = idx_b.query(ga.seq(na));
            for (NodeId nb : candidates) {
                float s = seq_similarity(ga.seq(na), gb.seq(nb));
                if (s >= params.min_sim)
                    cands.push_back({na, nb, s});
            }
        }

        // Greedy best-first assignment: sort by sim desc, assign each pair
        // if neither endpoint has been used yet.
        std::sort(cands.begin(), cands.end(),
                  [](const Candidate& a, const Candidate& b) {
                      return a.sim > b.sim;
                  });

        std::vector<bool> used_a(ga.num_nodes(), false);
        for (const auto& c : cands) {
            if (!used_a[c.na] && !used_b[c.nb]) {
                anchors.push_back({c.na, c.nb, c.sim});
                used_a[c.na] = used_b[c.nb] = true;
            }
        }
    }

    return anchors;
}

} // namespace g2g
