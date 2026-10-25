/**
 * src/graph.cpp
 * =============
 * VariationGraph builder methods - compact CSR layout.
 * TODO: implement add_node, add_path, finalise (edge-list → CSR).
 */
#include <g2g/graph.hpp>
#include <stdexcept>
#include <algorithm>

namespace g2g {

NodeId VariationGraph::add_node(const std::string& external_id,
                                 std::string_view sequence) {
    if (finalised_) throw std::logic_error("Cannot add nodes after finalise()");
    NodeId id = static_cast<NodeId>(node_seq_offsets.size());
    id_map[external_id] = id;
    node_names.push_back(external_id);
    node_seq_offsets.push_back(static_cast<uint64_t>(seq_data.size()));
    node_seq_lengths.push_back(static_cast<uint16_t>(sequence.size()));
    seq_data.insert(seq_data.end(), sequence.begin(), sequence.end());
    return id;
}

void VariationGraph::add_edge(NodeId u, NodeId v) {
    if (finalised_) throw std::logic_error("Cannot add edges after finalise()");
    edge_list_.emplace_back(u, v);
}

void VariationGraph::add_path(const std::string& pname,
                               const std::vector<NodeId>& nodes) {
    if (finalised_) throw std::logic_error("Cannot add paths after finalise()");
    path_names.push_back(pname);
    path_offsets.push_back(static_cast<uint32_t>(path_node_ids.size()));
    path_node_ids.insert(path_node_ids.end(), nodes.begin(), nodes.end());
    for (size_t i = 0; i + 1 < nodes.size(); ++i)
        add_edge(nodes[i], nodes[i+1]);
}

void VariationGraph::finalise() {
    if (finalised_) return;
    finalised_ = true;

    // Seal path offsets
    path_offsets.push_back(static_cast<uint32_t>(path_node_ids.size()));

    // Deduplicate edges: GFA files supply both explicit L lines AND implicit
    // edges from P/W lines (add_path calls add_edge for consecutive nodes).
    // Without deduplication, node degrees are inflated and the backbone
    // heuristic in find_anchors produces zero results.
    std::sort(edge_list_.begin(), edge_list_.end());
    edge_list_.erase(std::unique(edge_list_.begin(), edge_list_.end()),
                     edge_list_.end());

    uint32_t n = num_nodes();
    // Build out-edge CSR
    adj_out_ptr.assign(n + 1, 0);
    for (auto [u, v] : edge_list_) adj_out_ptr[u + 1]++;
    for (uint32_t i = 1; i <= n; ++i) adj_out_ptr[i] += adj_out_ptr[i-1];
    adj_out_data.resize(adj_out_ptr[n]);
    {
        std::vector<uint32_t> pos(adj_out_ptr.begin(), adj_out_ptr.begin() + n);
        for (auto [u, v] : edge_list_) adj_out_data[pos[u]++] = v;
    }
    // Build in-edge CSR
    adj_in_ptr.assign(n + 1, 0);
    for (auto [u, v] : edge_list_) adj_in_ptr[v + 1]++;
    for (uint32_t i = 1; i <= n; ++i) adj_in_ptr[i] += adj_in_ptr[i-1];
    adj_in_data.resize(adj_in_ptr[n]);
    {
        std::vector<uint32_t> pos(adj_in_ptr.begin(), adj_in_ptr.begin() + n);
        for (auto [u, v] : edge_list_) adj_in_data[pos[v]++] = u;
    }
    edge_list_.clear();
    edge_list_.shrink_to_fit();
}

} // namespace g2g
