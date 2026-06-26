#pragma once
/**
 * include/g2g/decompose.hpp
 * =========================
 * Snarl decomposition and snarl-tree layout.
 *
 * SnarlTree stores snarls in van Emde Boas (vEB) recursive order
 * (§2.5.2 of the paper) so that bottom-up DP traversal exhibits
 * cache-oblivious spatial locality.
 *
 * Two back-ends are provided at compile time:
 *   G2G_HAVE_HANDLEGRAPH  → delegate to libsnarls from vg toolkit.
 *                            Handles nested snarls, unary snarls,
 *                            cycles (via DAG projection).
 *   (default)             → built-in dominator approximation
 *                            (port of decompose.py), suitable for
 *                            DAGs only.  Used for testing.
 *
 * The vEB layout is applied by SnarlTree::build_veb_order() regardless
 * of which back-end produced the snarls.
 */

#include <g2g/types.hpp>
#include <g2g/graph.hpp>
#include <optional>   // add this line
#include <cstdint>
#include <utility>
#include <vector>
#include <span>

namespace g2g {

// ── Snarl tree ────────────────────────────────────────────────

/**
 * Flat array of Snarl objects in van Emde Boas order.
 * The tree structure is encoded in Snarl::children_begin/end —
 * all children of snarl[i] are stored at snarls[children_begin..children_end).
 *
 * vEB ordering guarantees that for any subtree of height h,
 * all its nodes fit in O(N^{1/2}) cache lines, giving an
 * O(N log_B N / B) cache complexity for bottom-up traversal
 * (optimal for comparison-based tree algorithms).
 */
/**
 * Per-snarl allele-path record in the flat path store.
 * Each snarl owns a contiguous block of allele_count() entries in
 * SnarlTree::allele_paths_data, indexed by
 *   [allele_offsets[i] .. allele_offsets[i] + allele_lengths[i])
 * where i is the flat index of the allele (0..allele_count()-1) within
 * this snarl's block starting at snarl.alleles_begin.
 * Each entry in allele_paths_data is itself a flat slice of NodeIds.
 */
struct SnarlTree {
    std::vector<Snarl>    snarls;     // vEB-ordered
    std::vector<uint32_t> roots;      // indices of top-level snarls
    uint32_t              depth{0};   // tree height

    // Flat allele-path storage (P3-6).
    // snarls[i].alleles_begin .. snarls[i].alleles_end index into
    // allele_path_offsets / allele_path_lengths; each allele's NodeIds
    // live at allele_paths_data[allele_path_offsets[k] .. +allele_path_lengths[k]].
    std::vector<uint32_t> allele_path_offsets;   // one entry per allele across all snarls
    std::vector<uint32_t> allele_path_lengths;   // one entry per allele
    std::vector<NodeId>   allele_paths_data;     // all NodeId sequences concatenated

    uint32_t size() const { return static_cast<uint32_t>(snarls.size()); }

    // Return the k-th allele path of snarl with alleles_begin=base.
    // allele_index is 0-based within the snarl (0 .. allele_count()-1).
    std::span<const NodeId> allele_path(uint32_t alleles_begin,
                                         uint32_t allele_index) const {
        uint32_t flat = alleles_begin + allele_index;
        return {allele_paths_data.data() + allele_path_offsets[flat],
                allele_path_lengths[flat]};
    }
};

// ── Decomposer ────────────────────────────────────────────────

class SnarlDecomposer {
public:
    explicit SnarlDecomposer(const VariationGraph& vg) : vg_(vg) {}

    /**
     * Decompose the graph into a SnarlTree in vEB order.
     * Raises std::runtime_error if the graph contains cycles
     * and G2G_HAVE_HANDLEGRAPH is not defined.
     */
    SnarlTree decompose() const;

private:
    const VariationGraph& vg_;

    // Build the flat SnarlTree from (source, sink) boundaries: enumerate
    // alleles via DFS, drop nested (non-top-level) snarls, annotate work, and
    // apply the vEB layout.  Shared by the built-in and vg boundary back-ends.
    static SnarlTree build_tree_from_boundaries(
        const VariationGraph& vg,
        const std::vector<std::pair<NodeId, NodeId>>& boundaries,
        const std::vector<uint32_t>& topo_pos);

    // Convert a list of snarls (DFS order) to vEB order in-place.
    static void apply_veb_order(SnarlTree& tree);

    // Compute subtree_work for each snarl (sum of k^3 over all descendants).
    static void annotate_work(SnarlTree& tree);
};

// ── External snarl back-end (snarls_vg.cpp) ───────────────────
//
// Ask the `vg` toolkit for snarl boundaries (source, sink) in NodeId space,
// shelling out to `vg snarls | vg view -Rj`.  Returns std::nullopt if vg is
// unavailable or yields nothing, so callers can fall back to the built-in
// decomposer.  Enabled at runtime via G2G_SNARLS=vg; vg located via $G2G_VG
// or PATH.  See snarls_vg.cpp.
std::optional<std::vector<std::pair<NodeId, NodeId>>>
vg_snarl_boundaries(const VariationGraph& vg);

// ── Topology helpers ──────────────────────────────────────────

/** Topological sort of a DAG; throws if cycles detected. */
std::vector<NodeId> topological_sort(const VariationGraph& vg);

/**
 * Cache the topological order in the graph for reuse.
 * Phase 3 fix: Python prototype called topological_sort 1.8M times;
 * C++ implementation computes it once per graph and stores it here.
 */
struct TopoCache {
    const VariationGraph*  vg{nullptr};
    std::vector<NodeId>    order;

    void build(const VariationGraph& g);
    const std::vector<NodeId>& get(const VariationGraph& g);
};

} // namespace g2g
