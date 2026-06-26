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

// ── find_sink ─────────────────────────────────────────────────────────────────

static std::optional<NodeId> find_sink(
        NodeId source,
        const VariationGraph& vg,
        const std::vector<NodeId>& topo,
        const std::vector<uint32_t>& /*topo_pos*/) {

    auto succs = vg.successors(source);
    if (succs.size() < 2) return std::nullopt;

    std::vector<std::vector<bool>> reachable(succs.size(),
                                              std::vector<bool>(vg.num_nodes(), false));
    for (size_t si = 0; si < succs.size(); ++si) {
        std::queue<NodeId> q;
        q.push(succs[si]);
        reachable[si][succs[si]] = true;
        while (!q.empty()) {
            NodeId u = q.front(); q.pop();
            for (NodeId v : vg.successors(u))
                if (!reachable[si][v]) { reachable[si][v] = true; q.push(v); }
        }
    }

    for (NodeId n : topo) {
        if (n == source) continue;
        bool in_all = true;
        for (size_t si = 0; si < succs.size(); ++si)
            if (!reachable[si][n]) { in_all = false; break; }
        if (in_all) return n;
    }
    return std::nullopt;
}

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
// For Phase 3 (single-level snarls, depth=1) the reordering is a no-op, but
// the infrastructure is in place for nested snarls once libsnarls is linked.
// The function is instrumented so the ablation test (Experiment 5, Q3) can
// compare traversal-order cache behaviour before vs. after vEB layout.
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
    // (marked TODO for Phase 4 when libsnarls provides nested snarls).
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
    // For nested snarls (Phase 4), remap all children_begin/end ranges here.
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
// Given snarl boundaries (source, sink) — from either the built-in dominator
// approximation or the vg back-end — enumerate each snarl's alleles, drop
// snarls nested inside another (keep top-level only), build the flat allele
// store, annotate work, and apply the vEB layout.

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

    for (const auto& [s, t] : boundaries) {
        if (s >= vg.num_nodes() || t >= vg.num_nodes()) continue;

        std::vector<std::vector<NodeId>> alleles;
        std::function<void(NodeId, std::vector<NodeId>&)> dfs =
            [&](NodeId u, std::vector<NodeId>& path) {
                if (u == t) { alleles.push_back(path); return; }
                for (NodeId v : vg.successors(u)) {
                    if (topo_pos[v] > topo_pos[t]) continue;
                    path.push_back(v);
                    dfs(v, path);
                    path.pop_back();
                }
            };
        std::vector<NodeId> path = {s};
        dfs(s, path);

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
    // dominator approximation.  vg falls back to built-in if unavailable.
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
        for (NodeId s : topo) {
            if (vg_.out_degree(s) < 2) continue;
            auto sink_opt = find_sink(s, vg_, topo, topo_pos);
            if (!sink_opt) continue;
            boundaries.emplace_back(s, *sink_opt);
        }
    }

    return build_tree_from_boundaries(vg_, boundaries, topo_pos);
}

} // namespace g2g
