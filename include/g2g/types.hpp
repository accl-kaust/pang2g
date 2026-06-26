#pragma once
/**
 * include/g2g/types.hpp
 * =====================
 * Core data structures. No computation here — mirrors types.py.
 *
 * Memory layout notes (for Phase 3 cache-oblivious work):
 *   - Node sequences are stored in VariationGraph::seq_data (flat byte array)
 *     and referenced by [offset, length] pairs.  No per-node std::string.
 *   - SnarlTree uses van Emde Boas layout (Phase 3, §2.5.2) — index helpers
 *     are in decompose.hpp.
 */

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace g2g {

// ── Node / Edge types ─────────────────────────────────────────

using NodeId = uint32_t;
constexpr NodeId kInvalidNode = ~NodeId{0};

struct NodeRef {
    uint32_t offset;   // byte offset into VariationGraph::seq_data
    uint16_t length;   // sequence length in bytes
    NodeId   id;       // stable external id
};

// ── Snarl ─────────────────────────────────────────────────────

/**
 * A snarl is a (source, sink) bubble.
 * allele_path_offsets: index into VariationGraph::path_data; each allele
 *   is a slice [begin, end) of NodeIds.
 * children_begin/end: range into the flat snarl array (vEB order).
 */
struct Snarl {
    NodeId   source;
    NodeId   sink;

    // Allele paths stored as flat slices into VariationGraph::path_node_ids
    uint32_t alleles_begin;   // first allele index
    uint32_t alleles_end;     // one-past-last allele index
    uint32_t allele_stride;   // max path length (paths zero-padded)

    // Children in snarl tree (indices into SnarlTree::snarls)
    uint32_t children_begin;
    uint32_t children_end;

    // Work estimate for scheduler (sum of k^3 over subtree)
    uint64_t subtree_work;

    uint32_t allele_count() const { return alleles_end - alleles_begin; }
    uint32_t child_count()  const { return children_end - children_begin; }
};

// ── Alignment records ─────────────────────────────────────────

struct AlleleMatch {
    uint16_t idx_a;
    uint16_t idx_b;
    float    sim;
};

struct SnarlAlignment {
    uint32_t snarl_a_idx;   // index into SnarlTree for graph A
    uint32_t snarl_b_idx;   // index into SnarlTree for graph B

    float    score;
    float    w_seq;
    float    w_topo;
    float    w_struct;

    std::vector<AlleleMatch> allele_matches;
    std::vector<uint16_t>    a_only;   // unmatched allele indices in snarl_a
    std::vector<uint16_t>    b_only;   // unmatched allele indices in snarl_b
};

struct AnchorMatch {
    NodeId node_a;
    NodeId node_b;
    float  sim;
};

// (snarl_source_name, allele_path) pairs unique to one graph
struct DeltaAllele {
    std::string          snarl_id;    // "<source>><sink>" string label
    std::vector<NodeId>  path;        // full allele path (source..sink inclusive)
    std::string          sequence;    // alternate sequence (ALT)

    // Reference coordinates (set by the reference-coordinate diff; -1 = unknown,
    // e.g. the snarl-alignment path which has no reference frame).
    std::string          chrom;       // reference path name
    int64_t              pos{-1};     // 0-based reference offset of the site
    std::string          ref;         // reference allele (REF); empty = insertion
};

struct AlignmentResult {
    std::vector<AnchorMatch>    anchors;
    std::vector<SnarlAlignment> snarl_alignments;
    float                       total_score{0.f};

    // Delta graph: alleles present in only one of the two graphs.
    // Populated by GraphAligner::build_deltas() (P3-7).
    std::vector<DeltaAllele>    a_only_alleles;   // in A, not in B
    std::vector<DeltaAllele>    b_only_alleles;   // in B, not in A
};

} // namespace g2g
