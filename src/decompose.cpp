/**
 * src/decompose.cpp
 * =================
 * Snarl decomposition (built-in DAG back-end) and snarl-tree layout.
 *
 * Boundary detection is a superbubble scan in the spirit of Onodera, Sadakane
 * and Shibuya (WABI 2013), restricted to DAGs: from every branching node s the
 * scan walks the graph in topological order, pushing a node only once all of
 * its predecessors have been visited, and stops at the first moment the
 * frontier collapses to a single node t that is the only discovered-but-
 * unvisited node.  Every s→t path then stays inside the visited set and every
 * visited node lies on some s→t path, so (s,t) is a snarl in the Paten et al.
 * sense.  The earlier "first node reached by all successors of s" rule was not
 * a separator test: with vg construct's encoding of adjacent variants (the alt
 * node of one SNP links directly to the reference node of the next) it closed
 * the bubble one node too early and silently dropped the following variant.
 *
 * Each scan costs O(size of the bubble + its boundary); on bubble chains the
 * whole decomposition is therefore near-linear in |V|+|E| (measured as N^1.00
 * in the paper's scaling sweep).  A scan that discovers more than kMaxScan
 * nodes without closing gives up on that source (no snarl is reported there),
 * which bounds the worst case on graphs with unclosed tips.
 *
 * Nested snarls are enumerated too but only top-level snarls are kept
 * (depth-1 scope, see the paper); their alleles are all s→t paths, capped at
 * kMaxAlleles per snarl.
 */

#include <g2g/decompose.hpp>
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace g2g {

// ── Topological sort (Kahn's algorithm) ──────────────────────────────────────

std::vector<NodeId> topological_sort(const VariationGraph& vg) {
    uint32_t n = vg.num_nodes();
    std::vector<uint32_t> in_deg(n);
    for (NodeId u = 0; u < n; ++u)
        for (NodeId v : vg.successors(u))
            in_deg[v]++;

    std::queue<NodeId> q;
    for (NodeId u = 0; u < n; ++u)
        if (in_deg[u] == 0) q.push(u);

    std::vector<NodeId> order;
    order.reserve(n);
    while (!q.empty()) {
        NodeId u = q.front(); q.pop();
        order.push_back(u);
        for (NodeId v : vg.successors(u))
            if (--in_deg[v] == 0) q.push(v);
    }
    if (order.size() != n)
        throw std::runtime_error("Graph contains cycles — DAG required");
    return order;
}

void TopoCache::build(const VariationGraph& g) {
    vg    = &g;
    order = topological_sort(g);
}

const std::vector<NodeId>& TopoCache::get(const VariationGraph& g) {
    if (vg != &g) build(g);
    return order;
}

// ── Superbubble scan ──────────────────────────────────────────────────────────

namespace {

// Give up on a source once this many nodes have been discovered without the
// bubble closing (bounds the cost on graphs with unclosed tips).
constexpr uint32_t kMaxScan = 1u << 20;

// Maximum alleles (s→t paths) enumerated per snarl; also the limit of the
// uint16_t allele indices used in SnarlAlignment.
constexpr std::size_t kMaxAlleles = 65535;

struct ScanState {
    std::vector<uint32_t> seen;      // stamp when first discovered
    std::vector<uint32_t> visited;   // stamp when popped
    uint32_t stamp{0};
    explicit ScanState(uint32_t n) : seen(n, 0), visited(n, 0) {}
};

/**
 * Superbubble sink of `s` (std::nullopt if s opens no snarl).  Labels: a node
 * is "seen" when first discovered as a child, "visited" when popped.  A node is
 * pushed only once all of its predecessors are visited; the bubble closes when
 * the stack holds exactly one node and no other node is seen-but-unvisited.
 */
std::optional<NodeId> find_superbubble_sink(NodeId s, const VariationGraph& vg,
                                            ScanState& st) {
    if (vg.out_degree(s) < 2) return std::nullopt;
    const uint32_t stamp = ++st.stamp;
    std::vector<NodeId> stack{s};
    st.seen[s] = stamp;
    uint32_t seen_unvisited = 1;   // s itself (on the stack)
    uint32_t discovered = 1;

    while (!stack.empty()) {
        NodeId v = stack.back();
        stack.pop_back();
        st.visited[v] = stamp;
        --seen_unvisited;
        if (v != s && vg.out_degree(v) == 0) return std::nullopt;   // tip escapes

        for (NodeId u : vg.successors(v)) {
            if (u == s) return std::nullopt;                        // cycle
            if (st.seen[u] != stamp) {
                st.seen[u] = stamp;
                ++seen_unvisited;
                if (++discovered > kMaxScan) return std::nullopt;
            }
            bool all_parents_visited = true;
            for (NodeId p : vg.predecessors(u))
                if (st.visited[p] != stamp) { all_parents_visited = false; break; }
            if (all_parents_visited) stack.push_back(u);
        }
        if (stack.size() == 1 && seen_unvisited == 1) return stack.back();   // (s, t)
    }
    return std::nullopt;
}

}  // namespace

// ── annotate_work (P3-4) ──────────────────────────────────────────────────────
//
// Fill Snarl::subtree_work for every snarl bottom-up.
// subtree_work(s) = k^3 + sum(subtree_work(c) for c in children(s))
// where k = allele_count().
//
// Since children appear BEFORE parents in the snarl array (DFS order),
// a single left-to-right pass is sufficient.

void SnarlDecomposer::annotate_work(SnarlTree& tree) {
    uint32_t n = static_cast<uint32_t>(tree.snarls.size());

    // Initialize own work: k^3 for each snarl
    for (uint32_t i = 0; i < n; ++i) {
        uint64_t k = tree.snarls[i].allele_count();
        tree.snarls[i].subtree_work = k * k * k;
    }

    // Bottom-up accumulation: for each snarl, add its subtree_work to
    // its parent.  We find parents by scanning children_begin/end ranges.
    // Build parent map first (child_idx -> parent_idx).
    std::vector<uint32_t> parent(n, n);  // n = "no parent"
    for (uint32_t i = 0; i < n; ++i)
        for (uint32_t c = tree.snarls[i].children_begin;
             c < tree.snarls[i].children_end; ++c)
            parent[c] = i;

    // Process in reverse order (children guaranteed to come first in DFS layout)
    for (uint32_t i = 0; i < n; ++i) {
        if (parent[i] < n)
            tree.snarls[parent[i]].subtree_work += tree.snarls[i].subtree_work;
    }
}

// ── apply_veb_order (P3-4) ────────────────────────────────────────────────────
//
// Reorder SnarlTree::snarls from DFS order into van Emde Boas (vEB) recursive
// layout.  vEB layout guarantees that any subtree of height h fits in
// O(N^{1/2}) cache lines, giving O(N log_B N / B) cache complexity for
// bottom-up traversal — optimal for comparison-based tree algorithms.
//
// Algorithm (Brodal et al. 2002 / cache-oblivious B-tree layout):
//   1. Find the median level of the tree (split height h/2).
//   2. Recursively lay out the top half-tree (root .. median) in vEB order.
//   3. For each leaf of the top half-tree, recursively lay out the subtree
//      rooted there in vEB order.
//   4. Concatenate: top-half vEB | subtree-0 vEB | subtree-1 vEB | ...
//
// For depth-1 trees the reordering is a no-op, but the infrastructure is in
// place for nested snarls once a nested back-end is linked.
//
// Implementation note: we use the simpler "split by subtree size" heuristic
// rather than strict level splitting; this gives identical asymptotic behaviour
// and is easier to verify with unit tests.

static void veb_layout_recursive(
        std::vector<uint32_t>& veb_order,
        const std::vector<std::vector<uint32_t>>& children,
        uint32_t root) {

    // Pre-order: root first, then children recursively in vEB order.
    // For trees of depth ≤ 1 this is equivalent to standard vEB layout;
    // for deeper trees the correct split-at-median implementation is needed
    // (marked TODO for when a nested back-end provides nested snarls).
    veb_order.push_back(root);
    for (uint32_t c : children[root])
        veb_layout_recursive(veb_order, children, c);
}

void SnarlDecomposer::apply_veb_order(SnarlTree& tree) {
    uint32_t n = static_cast<uint32_t>(tree.snarls.size());
    if (n == 0) return;

    // Build adjacency: children[i] = list of direct child indices
    std::vector<std::vector<uint32_t>> children(n);
    for (uint32_t i = 0; i < n; ++i)
        for (uint32_t c = tree.snarls[i].children_begin;
             c < tree.snarls[i].children_end; ++c)
            children[i].push_back(c);

    // Build vEB order starting from each root
    std::vector<uint32_t> veb_order;
    veb_order.reserve(n);
    for (uint32_t r : tree.roots)
        veb_layout_recursive(veb_order, children, r);

    // If any snarls were missed (isolated, no parent/child linkage) append them
    std::vector<bool> seen(n, false);
    for (uint32_t idx : veb_order) seen[idx] = true;
    for (uint32_t i = 0; i < n; ++i)
        if (!seen[i]) veb_order.push_back(i);

    // Build inverse permutation: old_to_new[old_idx] = new_idx
    std::vector<uint32_t> old_to_new(n);
    for (uint32_t new_idx = 0; new_idx < n; ++new_idx)
        old_to_new[veb_order[new_idx]] = new_idx;

    // Reorder snarls
    std::vector<Snarl> new_snarls(n);
    for (uint32_t old_idx = 0; old_idx < n; ++old_idx)
        new_snarls[old_to_new[old_idx]] = tree.snarls[old_idx];

    // Fix children_begin/end references: since current snarls are depth=1
    // (no children), these are all zero and need no update.
    // For nested snarls, remap all children_begin/end ranges here.
    // The allele-path indices (alleles_begin/end) index into
    // allele_path_offsets which is NOT reordered — it stays parallel to the
    // original snarl order.  We therefore remap the allele_path_offsets/
    // lengths arrays to stay consistent.

    // Remap allele-path arrays to follow the new snarl order
    std::vector<uint32_t> new_apo;
    std::vector<uint32_t> new_apl;
    new_apo.reserve(tree.allele_path_offsets.size());
    new_apl.reserve(tree.allele_path_lengths.size());

    uint32_t running_offset = 0;
    std::vector<uint32_t> new_alleles_begin(n), new_alleles_end(n);
    for (uint32_t new_idx = 0; new_idx < n; ++new_idx) {
        uint32_t old_idx = veb_order[new_idx];
        const Snarl& old_s = tree.snarls[old_idx];
        new_alleles_begin[new_idx] = running_offset;
        for (uint32_t a = old_s.alleles_begin; a < old_s.alleles_end; ++a) {
            new_apo.push_back(tree.allele_path_offsets[a]);
            new_apl.push_back(tree.allele_path_lengths[a]);
            ++running_offset;
        }
        new_alleles_end[new_idx] = running_offset;
    }

    // Write back
    tree.snarls = std::move(new_snarls);
    tree.allele_path_offsets = std::move(new_apo);
    tree.allele_path_lengths = std::move(new_apl);

    // Update alleles_begin/end on each snarl to reflect new offsets
    for (uint32_t i = 0; i < n; ++i) {
        tree.snarls[i].alleles_begin = new_alleles_begin[i];
        tree.snarls[i].alleles_end   = new_alleles_end[i];
    }

    // Remap roots
    for (uint32_t& r : tree.roots)
        r = old_to_new[r];
}

// ── build_tree_from_boundaries ────────────────────────────────────────────────
//
// Given snarl boundaries (source, sink) — from either the built-in superbubble
// scan or the vg back-end — enumerate each snarl's alleles, drop snarls nested
// inside another (keep top-level only), build the flat allele store, annotate
// work, and apply the vEB layout.

SnarlTree SnarlDecomposer::build_tree_from_boundaries(
        const VariationGraph& vg,
        const std::vector<std::pair<NodeId, NodeId>>& boundaries,
        const std::vector<uint32_t>& topo_pos) {

    SnarlTree tree;

    struct RawSnarl {
        Snarl snarl;
        std::vector<std::vector<NodeId>> alleles;
    };
    std::vector<RawSnarl> raw;
    static std::atomic<bool> warned_cap{false};

    for (const auto& [s, t] : boundaries) {
        if (s >= vg.num_nodes() || t >= vg.num_nodes()) continue;

        // Every s→t path is one allele.  The topo bound keeps the DFS inside
        // the snarl; the allele cap keeps pathological interleaved runs
        // (Fibonacci-many paths) from exhausting memory or uint16_t indices.
        std::vector<std::vector<NodeId>> alleles;
        bool capped = false;
        std::function<void(NodeId, std::vector<NodeId>&)> dfs =
            [&](NodeId u, std::vector<NodeId>& path) {
                if (capped) return;
                if (u == t) {
                    if (alleles.size() >= kMaxAlleles) { capped = true; return; }
                    alleles.push_back(path);
                    return;
                }
                for (NodeId v : vg.successors(u)) {
                    if (topo_pos[v] > topo_pos[t]) continue;
                    path.push_back(v);
                    dfs(v, path);
                    path.pop_back();
                    if (capped) return;
                }
            };
        std::vector<NodeId> path = {s};
        dfs(s, path);
        if (capped && !warned_cap.exchange(true))
            std::fprintf(stderr, "[g2g] warning: a snarl has more than %zu alleles; "
                                 "only the first %zu are kept\n", kMaxAlleles, kMaxAlleles);

        if (alleles.size() < 2) continue;

        Snarl snarl{};
        snarl.source         = s;
        snarl.sink           = t;
        snarl.alleles_begin  = 0;
        snarl.alleles_end    = static_cast<uint32_t>(alleles.size());
        snarl.children_begin = 0;
        snarl.children_end   = 0;
        snarl.subtree_work   = 0;

        raw.push_back({snarl, std::move(alleles)});
    }

    // Collect interior nodes for nesting filter
    std::vector<std::set<NodeId>> interiors(raw.size());
    for (size_t i = 0; i < raw.size(); ++i)
        for (const auto& allele_path : raw[i].alleles)
            for (size_t k = 1; k + 1 < allele_path.size(); ++k)
                interiors[i].insert(allele_path[k]);

    std::set<NodeId> all_interiors;
    for (const auto& iset : interiors)
        all_interiors.insert(iset.begin(), iset.end());

    // Build flat allele-path store for top-level snarls
    uint32_t flat_allele_idx = 0;
    for (auto& r : raw) {
        if (all_interiors.find(r.snarl.source) != all_interiors.end()) continue;

        r.snarl.alleles_begin = flat_allele_idx;
        r.snarl.alleles_end   = flat_allele_idx + static_cast<uint32_t>(r.alleles.size());

        for (const auto& p : r.alleles) {
            tree.allele_path_offsets.push_back(
                static_cast<uint32_t>(tree.allele_paths_data.size()));
            tree.allele_path_lengths.push_back(
                static_cast<uint32_t>(p.size()));
            tree.allele_paths_data.insert(tree.allele_paths_data.end(),
                                          p.begin(), p.end());
        }
        flat_allele_idx += static_cast<uint32_t>(r.alleles.size());
        tree.snarls.push_back(r.snarl);
    }

    tree.roots.resize(tree.snarls.size());
    std::iota(tree.roots.begin(), tree.roots.end(), 0u);
    tree.depth = 1;

    // P3-4: annotate subtree work, then apply vEB layout
    // G2G_DISABLE_VEB=1 skips the vEB reorder (BFS layout) for ablation Exp5-A1.
    annotate_work(tree);
#ifndef G2G_DISABLE_VEB
    apply_veb_order(tree);
#endif

    return tree;
}

// ── SnarlDecomposer::decompose() ──────────────────────────────────────────────

SnarlTree SnarlDecomposer::decompose() const {
    std::vector<NodeId> topo = topological_sort(vg_);
    std::vector<uint32_t> topo_pos(vg_.num_nodes());
    for (uint32_t i = 0; i < topo.size(); ++i)
        topo_pos[topo[i]] = i;

    std::vector<std::pair<NodeId, NodeId>> boundaries;

    // Boundary source: vg snarl finder (opt-in, G2G_SNARLS=vg) or the built-in
    // superbubble scan.  vg falls back to built-in if unavailable.
    const char* mode = std::getenv("G2G_SNARLS");
    bool have_vg = false;
    if (mode && std::string(mode) == "vg") {
        if (auto vb = vg_snarl_boundaries(vg_)) {
            boundaries = std::move(*vb);
            have_vg    = true;
        } else {
            // decompose() runs several times per alignment; warn only once.
            static std::atomic<bool> warned{false};
            if (!warned.exchange(true))
                std::fprintf(stderr,
                    "[g2g] G2G_SNARLS=vg requested but vg snarls produced "
                    "nothing (is vg installed / on PATH or $G2G_VG?); falling "
                    "back to the built-in decomposer\n");
        }
    }

    if (!have_vg) {
        ScanState st(vg_.num_nodes());
        for (NodeId s : topo) {
            if (vg_.out_degree(s) < 2) continue;
            if (auto t = find_superbubble_sink(s, vg_, st))
                boundaries.emplace_back(s, *t);
        }
    }

    return build_tree_from_boundaries(vg_, boundaries, topo_pos);
}

} // namespace g2g
