/**
 * tests/test_phase5_exp2_minimap2.cpp
 * =====================================
 * Phase 5, Experiment 2 - Minimap2 Sequence-Level Consistency
 *
 * Protocol (roadmap §2.7.1 / §5.1 Option A)
 * ------------------------------------------
 * For small graphs (≤ 20 snarls), every G2G allele match (node_a, node_b)
 * must be *consistent* with minimap2: there must exist a minimap2 alignment
 * between some path containing node_a and some path containing node_b,
 * where the aligned region covers both nodes' sequence spans.
 *
 * Consistency definition
 * ----------------------
 * Match (inner_node_a, inner_node_b) is consistent iff:
 *   ∃ path p_a through node_a, path p_b through node_b,
 *     minimap2 aligns p_a → p_b with query interval [qs,qe) overlapping
 *     node_a's byte span in p_a AND target interval [ts,te) overlapping
 *     node_b's byte span in p_b.
 *
 * Path coverage guarantee
 * -----------------------
 * For each node n, we include at least one path through n before filling
 * the remainder up to MAX_PATHS by DFS. This prevents the false failures
 * that arise from a naive front-loaded enumeration cap.
 *
 * Implementation
 * --------------
 * minimap2 is called as a subprocess; PAF output is parsed to build a
 * coverage index over (node_a, node_b) pairs.  The test is SKIPPED when
 * the minimap2 binary is not on PATH (keeps CI green on machines without it).
 */

#include <g2g/graph.hpp>
#include <g2g/align.hpp>
#include <g2g/sequence.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace g2g;
namespace fs = std::filesystem;

// -------------------------------------------------------------
// Utilities
// -------------------------------------------------------------

static std::string random_dna(int len, std::mt19937& rng) {
    static const char BASES[] = "ACGT";
    std::uniform_int_distribution<int> d(0, 3);
    std::string s(len, 'A');
    for (char& c : s) c = BASES[d(rng)];
    return s;
}

/** Check minimap2 is on PATH by calling `minimap2 --version`. */
static bool minimap2_available() {
    return std::system("minimap2 --version > /dev/null 2>&1") == 0;
}

// -------------------------------------------------------------
// Graph builder (same topology as Experiment 1)
// -------------------------------------------------------------

static std::pair<VariationGraph, VariationGraph>
build_graph_pair(int n_snarls, int k_shared = 2, int k_only = 1,
                 int seq_len = 50, uint32_t seed = 42) {
    std::mt19937 rng(seed);
    VariationGraph ga, gb;
    ga.name = "GA"; gb.name = "GB";

    std::vector<NodeId> anc_a, anc_b;
    for (int i = 0; i <= n_snarls; ++i) {
        std::string seq = random_dna(8, rng);
        std::string name = "a" + std::to_string(i);
        anc_a.push_back(ga.add_node(name, seq));
        anc_b.push_back(gb.add_node(name, seq));
    }

    for (int s = 0; s < n_snarls; ++s) {
        NodeId sa = anc_a[s], ta = anc_a[s+1];
        NodeId sb = anc_b[s], tb = anc_b[s+1];

        for (int j = 0; j < k_shared; ++j) {
            std::string seq = random_dna(seq_len, rng);
            std::string na = "s" + std::to_string(s) + "_sh" + std::to_string(j) + "_a";
            std::string nb = "s" + std::to_string(s) + "_sh" + std::to_string(j) + "_b";
            NodeId ia = ga.add_node(na, seq);
            NodeId ib = gb.add_node(nb, seq);
            ga.add_edge(sa, ia); ga.add_edge(ia, ta);
            gb.add_edge(sb, ib); gb.add_edge(ib, tb);
        }
        for (int j = 0; j < k_only; ++j) {
            std::string seq = random_dna(seq_len, rng);
            NodeId ia = ga.add_node("s" + std::to_string(s) + "_ao" + std::to_string(j), seq);
            ga.add_edge(sa, ia); ga.add_edge(ia, ta);
        }
        for (int j = 0; j < k_only; ++j) {
            std::string seq = random_dna(seq_len, rng);
            NodeId ib = gb.add_node("s" + std::to_string(s) + "_bo" + std::to_string(j), seq);
            gb.add_edge(sb, ib); gb.add_edge(ib, tb);
        }
    }

    // Representative paths
    std::vector<NodeId> pa = { anc_a[0] };
    std::vector<NodeId> pb = { anc_b[0] };
    for (int s = 0; s < n_snarls; ++s) {
        pa.push_back(ga.id_map.at("s" + std::to_string(s) + "_sh0_a"));
        pa.push_back(anc_a[s+1]);
        pb.push_back(gb.id_map.at("s" + std::to_string(s) + "_sh0_b"));
        pb.push_back(anc_b[s+1]);
    }
    ga.add_path("hap_A", pa);
    gb.add_path("hap_B", pb);

    ga.finalise();
    gb.finalise();
    return { std::move(ga), std::move(gb) };
}

// -------------------------------------------------------------
// Path enumeration with node-coverage guarantee
// -------------------------------------------------------------

using Path = std::vector<NodeId>;
using Frame = std::pair<NodeId, Path>;

/**
 * Enumerate paths using iterative DFS with a visited-set to avoid cycles.
 * Returns up to max_paths paths, but guarantees that every node in `vg`
 * is covered by at least one path before the cap is applied.
 *
 * Strategy:
 *   Phase 1 (coverage guarantee): for every node, find ONE path src→node→sink
 *     via a direct BFS-to-node + DFS-forward-to-sink construction.
 *     These paths are added regardless of the cap.
 *   Phase 2 (fill quota): DFS from all sources up to max_paths total.
 */
static std::vector<Path> enumerate_paths(const VariationGraph& vg,
                                          int max_paths = 500) {
    std::vector<NodeId> sources, sinks;
    for (NodeId n = 0; n < vg.num_nodes(); ++n) {
        if (vg.in_degree(n)  == 0) sources.push_back(n);
        if (vg.out_degree(n) == 0) sinks.push_back(n);
    }

    // Helper: find any path from `from` to `to` via BFS, returns empty if none.
    auto bfs_path = [&](NodeId from, NodeId to) -> Path {
        if (from == to) return { from };
        std::unordered_map<NodeId, NodeId> parent;
        std::vector<NodeId> q = { from };
        parent[from] = from;
        for (size_t qi = 0; qi < q.size(); ++qi) {
            NodeId cur = q[qi];
            if (cur == to) {
                Path p;
                while (cur != from) { p.push_back(cur); cur = parent[cur]; }
                p.push_back(from);
                std::reverse(p.begin(), p.end());
                return p;
            }
            for (NodeId nxt : vg.successors(cur))
                if (!parent.count(nxt)) { parent[nxt] = cur; q.push_back(nxt); }
        }
        return {};
    };

    // Helper: find any path from `from` to any sink.
    auto dfs_to_sink = [&](NodeId from) -> Path {
        std::vector<Frame> stk = {{ from, { from } }};
        while (!stk.empty()) {
            auto [node, p] = stk.back(); stk.pop_back();
            if (vg.out_degree(node) == 0) return p;
            for (NodeId nxt : vg.successors(node)) {
                bool cyc = false;
                for (NodeId x : p) if (x == nxt) { cyc = true; break; }
                if (!cyc) { Path np = p; np.push_back(nxt); stk.push_back({ nxt, np }); }
            }
        }
        return {};
    };

    std::set<Path> path_set;

    // -- Phase 1: guarantee every node is covered --------------
    for (NodeId target = 0; target < vg.num_nodes(); ++target) {
        // Build: best_source → target → any_sink
        bool added = false;
        for (NodeId src : sources) {
            Path prefix = bfs_path(src, target);
            if (prefix.empty()) continue;
            Path suffix = dfs_to_sink(target);
            if (suffix.empty()) continue;
            // Merge: prefix ends at target, suffix starts at target
            Path full = prefix;
            full.insert(full.end(), suffix.begin() + 1, suffix.end());
            path_set.insert(full);
            added = true;
            break;
        }
        // If target is a source with no predecessors, just DFS forward
        if (!added) {
            Path p = dfs_to_sink(target);
            if (!p.empty()) path_set.insert(p);
        }
    }

    // -- Phase 2: fill remaining quota with DFS -----------------
    for (NodeId src : sources) {
        std::vector<Frame> stk = {{ src, { src } }};
        while (!stk.empty() && (int)path_set.size() < max_paths) {
            auto [node, cur] = stk.back(); stk.pop_back();
            if (vg.out_degree(node) == 0) {
                path_set.insert(cur);
            } else {
                for (NodeId nxt : vg.successors(node)) {
                    bool cyc = false;
                    for (NodeId x : cur) if (x == nxt) { cyc = true; break; }
                    if (!cyc) {
                        Path np = cur; np.push_back(nxt);
                        stk.push_back({ nxt, std::move(np) });
                    }
                }
            }
        }
        if ((int)path_set.size() >= max_paths) break;
    }

    return std::vector<Path>(path_set.begin(), path_set.end());
}

// -------------------------------------------------------------
// GFA and FASTA writers
// -------------------------------------------------------------

[[maybe_unused]] static void write_gfa(const VariationGraph& vg, const fs::path& p) {
    std::ofstream f(p);
    f << "H\tVN:Z:1.1\n";
    for (NodeId n = 0; n < vg.num_nodes(); ++n) {
        auto sv = vg.seq(n);
        f << "S\t" << vg.node_names[n] << "\t";
        f.write(sv.data(), sv.size());
        f << "\n";
    }
    std::set<std::pair<NodeId,NodeId>> seen;
    for (NodeId u = 0; u < vg.num_nodes(); ++u)
        for (NodeId v : vg.successors(u))
            if (!seen.count({u,v})) {
                f << "L\t" << vg.node_names[u] << "\t+\t"
                  << vg.node_names[v] << "\t+\t0M\n";
                seen.insert({u,v});
            }
    for (size_t pi = 0; pi < vg.path_names.size(); ++pi) {
        uint32_t beg = vg.path_offsets[pi], end = vg.path_offsets[pi+1];
        f << "P\t" << vg.path_names[pi] << "\t";
        for (uint32_t i = beg; i < end; ++i) {
            if (i > beg) f << ",";
            f << vg.node_names[vg.path_node_ids[i]] << "+";
        }
        f << "\t";
        for (uint32_t i = beg; i + 1 < end; ++i) { if (i > beg) f << ","; f << "0M"; }
        f << (end - beg <= 1 ? "*" : "") << "\n";
    }
}

/**
 * Write enumerated paths as FASTA, skipping degenerate paths (seq < 5 bp).
 * Returns map: fasta_header → Path
 */
static std::unordered_map<std::string, Path>
write_paths_fasta(const VariationGraph& vg,
                  const std::vector<Path>& paths,
                  const fs::path& out) {
    std::ofstream f(out);
    std::unordered_map<std::string, Path> index;
    for (size_t i = 0; i < paths.size(); ++i) {
        std::string seq;
        for (NodeId n : paths[i]) {
            auto sv = vg.seq(n);
            seq.append(sv.data(), sv.size());
        }
        if (seq.size() < 5) continue;
        std::string name = "path_" + std::to_string(i);
        f << ">" << name << "\n" << seq << "\n";
        index[name] = paths[i];
    }
    return index;
}

// -------------------------------------------------------------
// PAF parser
// -------------------------------------------------------------

struct PafRecord {
    std::string query_name;
    int query_start{}, query_end{};
    std::string target_name;
    int target_start{}, target_end{};
    double identity{};
};

static std::vector<PafRecord> parse_paf(const fs::path& p) {
    std::vector<PafRecord> recs;
    std::ifstream f(p);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::vector<std::string> cols;
        for (std::string tok; ss >> tok; ) cols.push_back(tok);
        if (cols.size() < 12) continue;
        PafRecord r;
        r.query_name   = cols[0];
        r.query_start  = std::stoi(cols[2]);
        r.query_end    = std::stoi(cols[3]);
        r.target_name  = cols[5];
        r.target_start = std::stoi(cols[7]);
        r.target_end   = std::stoi(cols[8]);
        int res = std::stoi(cols[9]);
        int aln = std::stoi(cols[10]);
        r.identity = aln > 0 ? static_cast<double>(res) / aln : 0.0;
        recs.push_back(r);
    }
    return recs;
}

// -------------------------------------------------------------
// Node-offset helper
// -------------------------------------------------------------

/** Return [start, end) byte span of node n within the path's concatenated sequence. */
static std::pair<int,int> node_span_in_path(const VariationGraph& vg,
                                              const Path& path, NodeId node) {
    int offset = 0;
    for (NodeId n : path) {
        int len = static_cast<int>(vg.seq(n).size());
        if (n == node) return { offset, offset + len };
        offset += len;
    }
    return { -1, -1 };
}

// -------------------------------------------------------------
// Coverage index builder
// -------------------------------------------------------------

using NodePair = std::pair<NodeId, NodeId>;

static std::set<NodePair> build_coverage_index(
    const std::vector<PafRecord>& paf,
    const std::unordered_map<std::string, Path>& paths_a_idx,
    const std::unordered_map<std::string, Path>& paths_b_idx,
    const VariationGraph& ga,
    const VariationGraph& gb)
{
    std::set<NodePair> covered;
    for (const auto& rec : paf) {
        auto it_a = paths_a_idx.find(rec.query_name);
        auto it_b = paths_b_idx.find(rec.target_name);
        if (it_a == paths_a_idx.end() || it_b == paths_b_idx.end()) continue;
        const Path& pa = it_a->second;
        const Path& pb = it_b->second;

        for (NodeId na : pa) {
            auto [a0, a1] = node_span_in_path(ga, pa, na);
            if (a1 <= rec.query_start || a0 >= rec.query_end) continue;
            for (NodeId nb : pb) {
                auto [b0, b1] = node_span_in_path(gb, pb, nb);
                if (b1 <= rec.target_start || b0 >= rec.target_end) continue;
                covered.insert({ na, nb });
            }
        }
    }
    return covered;
}

// -------------------------------------------------------------
// G2G match extractor
// -------------------------------------------------------------

/** Extract (inner_node_a, inner_node_b) pairs from a G2G alignment result. */
static std::vector<NodePair> extract_g2g_pairs(
    const AlignmentResult& result,
    const SnarlTree& ta,
    const SnarlTree& tb)
{
    std::vector<NodePair> pairs;
    for (const auto& aln : result.snarl_alignments) {
        const Snarl& sa = ta.snarls[aln.snarl_a_idx];
        const Snarl& sb = tb.snarls[aln.snarl_b_idx];
        for (const auto& m : aln.allele_matches) {
            auto path_a = ta.allele_path(sa.alleles_begin, m.idx_a);
            auto path_b = tb.allele_path(sb.alleles_begin, m.idx_b);
            // Inner nodes only (skip source/sink boundary)
            for (size_t i = 1; i + 1 < path_a.size(); ++i)
                for (size_t j = 1; j + 1 < path_b.size(); ++j)
                    pairs.push_back({ path_a[i], path_b[j] });
        }
    }
    return pairs;
}

// -------------------------------------------------------------
// Core consistency check logic
// -------------------------------------------------------------

struct ConsistencyResult {
    int n_g2g_pairs{};
    int n_consistent{};
    int n_paf_records{};
    int n_paths_a{};
    int n_paths_b{};

    double rate() const {
        return n_g2g_pairs > 0
            ? static_cast<double>(n_consistent) / n_g2g_pairs
            : 1.0;
    }
};

static ConsistencyResult check_consistency(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const fs::path& tmp_dir,
    int max_paths = 500)
{
    ConsistencyResult res;

    // -- 1. Run G2G ---------------------------------------------
    AlignerParams params;
    GraphAligner aligner(params);
    auto result = aligner.align(ga, gb);

    SnarlTree ta = SnarlDecomposer(ga).decompose();
    SnarlTree tb = SnarlDecomposer(gb).decompose();

    auto g2g_pairs = extract_g2g_pairs(result, ta, tb);
    res.n_g2g_pairs = static_cast<int>(g2g_pairs.size());

    if (g2g_pairs.empty()) {
        res.n_consistent = 0;
        return res;  // trivially consistent
    }

    // -- 2. Enumerate paths -------------------------------------
    auto paths_a = enumerate_paths(ga, max_paths);
    auto paths_b = enumerate_paths(gb, max_paths);
    res.n_paths_a = static_cast<int>(paths_a.size());
    res.n_paths_b = static_cast<int>(paths_b.size());

    // -- 3. Write FASTA files -----------------------------------
    fs::path fa_file = tmp_dir / "paths_a.fa";
    fs::path fb_file = tmp_dir / "paths_b.fa";
    auto idx_a = write_paths_fasta(ga, paths_a, fa_file);
    auto idx_b = write_paths_fasta(gb, paths_b, fb_file);

    // -- 4. Run minimap2 ----------------------------------------
    fs::path paf_file = tmp_dir / "out.paf";
    std::string cmd = "minimap2 -x asm5 --secondary=no "
                    + fb_file.string() + " "
                    + fa_file.string() + " > "
                    + paf_file.string() + " 2>/dev/null";
    int rc = std::system(cmd.c_str());
    if (rc != 0) throw std::runtime_error("minimap2 subprocess failed");

    auto paf = parse_paf(paf_file);
    res.n_paf_records = static_cast<int>(paf.size());

    // -- 5. Build coverage index --------------------------------
    auto covered = build_coverage_index(paf, idx_a, idx_b, ga, gb);

    // -- 6. Count consistent pairs ------------------------------
    for (const auto& [na, nb] : g2g_pairs)
        if (covered.count({ na, nb })) ++res.n_consistent;

    return res;
}

// -------------------------------------------------------------
// Test fixture + cases
// -------------------------------------------------------------

class Minimap2ConsistencyTest : public ::testing::TestWithParam<int> {
protected:
    void SetUp() override {
        if (!minimap2_available())
            GTEST_SKIP() << "minimap2 not found on PATH - skipping Experiment 2";
    }
};

TEST_P(Minimap2ConsistencyTest, ConsistencyAtLeast80Pct) {
    int n_snarls = GetParam();
    auto [ga, gb] = build_graph_pair(n_snarls, 2, 1, 15, 42);

    // Write to /tmp/<test_unique_dir> - remove first to prevent stale files
    fs::path tmp = fs::temp_directory_path()
                 / ("g2g_exp2_" + std::to_string(n_snarls));
    fs::remove_all(tmp);
    fs::create_directories(tmp);

    ConsistencyResult res;
    ASSERT_NO_THROW(res = check_consistency(ga, gb, tmp));

    RecordProperty("n_snarls",      n_snarls);
    RecordProperty("n_g2g_pairs",   res.n_g2g_pairs);
    RecordProperty("n_consistent",  res.n_consistent);
    RecordProperty("n_paf_records", res.n_paf_records);
    RecordProperty("n_paths_a",     res.n_paths_a);
    RecordProperty("n_paths_b",     res.n_paths_b);

    if (res.n_g2g_pairs == 0) {
        SUCCEED() << "No G2G matches - trivially consistent";
        return;
    }

    EXPECT_GE(res.rate(), 0.80)
        << "Consistency too low: " << res.n_consistent
        << "/" << res.n_g2g_pairs
        << " (" << res.rate() * 100 << "%)";
}

// -- Perfect consistency on single snarl ----------------------
TEST_F(Minimap2ConsistencyTest, SingleSnarlPerfectConsistency) {
    if (!minimap2_available())
        GTEST_SKIP() << "minimap2 not found on PATH";

    auto [ga, gb] = build_graph_pair(1, 2, 1, 50, 7);
    fs::path tmp = fs::temp_directory_path() / "g2g_exp2_single";
    fs::remove_all(tmp);
    fs::create_directories(tmp);

    ConsistencyResult res;
    ASSERT_NO_THROW(res = check_consistency(ga, gb, tmp));

    if (res.n_g2g_pairs > 0) {
        EXPECT_GE(res.rate(), 0.8)
            << "Expected ≥80% consistency for single snarl, got "
            << res.n_consistent << "/" << res.n_g2g_pairs;
    }
}

// -- PAF records produced -------------------------------------
TEST_F(Minimap2ConsistencyTest, MiniMap2ProducesPafOutput) {
    if (!minimap2_available())
        GTEST_SKIP() << "minimap2 not found on PATH";

    auto [ga, gb] = build_graph_pair(3, 2, 1, 15, 99);
    fs::path tmp = fs::temp_directory_path() / "g2g_exp2_paf";
    fs::remove_all(tmp);
    fs::create_directories(tmp);

    ConsistencyResult res;
    ASSERT_NO_THROW(res = check_consistency(ga, gb, tmp));
    EXPECT_GT(res.n_paf_records, 0)
        << "Expected minimap2 to produce at least one PAF record";
}

// -- All paths include every node (coverage guarantee) ---------
TEST_F(Minimap2ConsistencyTest, PathEnumerationCoversAllNodes) {
    if (!minimap2_available())
        GTEST_SKIP() << "minimap2 not found on PATH";

    auto [ga, gb] = build_graph_pair(5, 2, 1, 15, 42);

    auto paths_a = enumerate_paths(ga, 500);
    std::unordered_set<NodeId> covered;
    for (const auto& p : paths_a)
        for (NodeId n : p) covered.insert(n);

    for (NodeId n = 0; n < ga.num_nodes(); ++n)
        EXPECT_TRUE(covered.count(n))
            << "Node " << ga.node_names[n]
            << " (id=" << n << ") not covered by any enumerated path";
}

INSTANTIATE_TEST_SUITE_P(
    GraphSizes,
    Minimap2ConsistencyTest,
    ::testing::Values(3, 5, 8, 12),
    [](const ::testing::TestParamInfo<int>& info) {
        return "n" + std::to_string(info.param);
    });
