#pragma once
/**
 * include/g2g/graph.hpp
 * =====================
 * VariationGraph: compact adjacency list with flat sequence storage.
 *
 * Memory layout (cache-friendly, avoids pointer chasing):
 *   seq_data          : flat byte array of all node sequences, concatenated
 *   node_seq_offsets  : [i] = byte offset of node i's sequence in seq_data
 *   node_seq_lengths  : [i] = byte length of node i's sequence
 *   adj_out           : CSR adjacency (out-edges): adj_targets[adj_out[i]..adj_out[i+1]]
 *   adj_in            : CSR adjacency (in-edges)
 *   path_node_ids     : flat array of node IDs for all haplotype paths
 *   path_offsets      : [i] = start index of path i in path_node_ids
 *
 * This is the data layout the C++ DP kernel will read; keeping sequences
 * in a single flat array maximises L2/L3 cache utilization during alignment.
 */

#include <g2g/types.hpp>

#include <cassert>
#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <g2g/span.hpp>   // portable g2g::span (falls back on GCC 9)
#include <initializer_list>

namespace g2g {

class VariationGraph {
public:
    std::string name;

    // -- Sequence storage --------------------------------------
    // All node sequences packed into one flat array.
    std::vector<char>     seq_data;
    std::vector<uint64_t> node_seq_offsets;  // size = num_nodes; 64-bit so a
                                             // single graph's seq_data may
                                             // exceed 4 GiB (multi-Gbp graphs)
    std::vector<uint16_t> node_seq_lengths;  // size = num_nodes

    // -- Adjacency (CSR format) ---------------------------------
    std::vector<uint32_t> adj_out_ptr;    // size = num_nodes + 1
    std::vector<NodeId>   adj_out_data;
    std::vector<uint32_t> adj_in_ptr;     // size = num_nodes + 1
    std::vector<NodeId>   adj_in_data;

    // -- Named paths -------------------------------------------
    std::vector<uint32_t>      path_offsets;   // size = num_paths + 1
    std::vector<NodeId>        path_node_ids;  // flat
    std::vector<std::string>   path_names;

    // -- External-id ↔ internal-index mapping -----------------
    // Python prototype uses string node IDs; we keep a map for GFA import.
    std::unordered_map<std::string, NodeId> id_map;
    std::vector<std::string>               node_names;  // reverse map

    // -- Accessors ---------------------------------------------
    uint32_t num_nodes() const { return static_cast<uint32_t>(node_seq_offsets.size()); }
    uint32_t num_paths() const { return static_cast<uint32_t>(path_names.size()); }

    std::string_view seq(NodeId n) const {
        assert(n < num_nodes());
        return {seq_data.data() + node_seq_offsets[n], node_seq_lengths[n]};
    }

    // Out-neighbours of node n
    g2g::span<const NodeId> successors(NodeId n) const {
        return {adj_out_data.data() + adj_out_ptr[n],
                adj_out_ptr[n + 1] - adj_out_ptr[n]};
    }
    // In-neighbours of node n
    g2g::span<const NodeId> predecessors(NodeId n) const {
        return {adj_in_data.data() + adj_in_ptr[n],
                adj_in_ptr[n + 1] - adj_in_ptr[n]};
    }

    uint32_t out_degree(NodeId n) const { return adj_out_ptr[n+1] - adj_out_ptr[n]; }
    uint32_t in_degree (NodeId n) const { return adj_in_ptr[n+1]  - adj_in_ptr[n]; }

    // Path node IDs as a span
    g2g::span<const NodeId> path_nodes(uint32_t path_idx) const {
        return {path_node_ids.data() + path_offsets[path_idx],
                path_offsets[path_idx + 1] - path_offsets[path_idx]};
    }

    // -- Builders (used by GFA reader and tests) ---------------
    // These produce the Python prototype's builder interface on top
    // of the compact CSR layout.  Call finalise() once after all
    // add_node/add_path calls.
    // -- Builders (used by GFA reader and tests) ---------------
    NodeId add_node(const std::string& external_id, std::string_view sequence);
    void   add_edge(NodeId u, NodeId v);
    void   add_path(const std::string& name, const std::vector<NodeId>& nodes);
    void   add_path(const std::string& name, std::initializer_list<NodeId> nodes) {
        add_path(name, std::vector<NodeId>(nodes));
    }
    void   finalise();

private:
    // Temporary edge list - converted to CSR by finalise()
    std::vector<std::pair<NodeId,NodeId>> edge_list_;
    bool finalised_{false};
};

} // namespace g2g
