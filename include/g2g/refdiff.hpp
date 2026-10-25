#pragma once
/**
 * include/g2g/refdiff.hpp
 * =======================
 * Reference-coordinate graph-to-graph diff.
 *
 * The snarl-alignment correspondence (correspond_snarls) pairs snarls by their
 * boundary node sequences.  That is brittle on real pangenome graphs: two
 * graphs built from the same reference but different variant sets segment the
 * shared reference into different nodes, so even shared variants have different
 * boundary sequences and fail to line up.
 *
 * reference_diff() instead anchors both graphs to a shared embedded reference
 * path and matches snarls by the reference position at which their bubble opens
 * (source_offset + source_length).  That position is a property of the
 * reference coordinate, independent of how each graph chops its nodes, so
 * corresponding variants line up regardless of segmentation.  At each position
 * the alternate alleles are compared by sequence: alleles in both graphs are
 * shared, alleles in only one go to that graph's delta - with real reference
 * coordinates (chrom/pos) attached.
 *
 * Requires both graphs to embed the reference as a path: a GFA P-line named
 * `ref_path`, or a W-line whose PanSN name (sample#hap#contig[:range]) has
 * `ref_path` as its sample or contig, e.g. --ref GRCh38 or --ref chr22.  A
 * reference spanning several contigs is diffed contig by contig.  Returns
 * std::nullopt if either graph lacks such a path.
 */

#include <g2g/types.hpp>
#include <g2g/graph.hpp>
#include <g2g/decompose.hpp>

#include <optional>
#include <string>

namespace g2g {

std::optional<AlignmentResult> reference_diff(
    const VariationGraph& ga,
    const VariationGraph& gb,
    const SnarlTree&      ta,
    const SnarlTree&      tb,
    const std::string&    ref_path);

} // namespace g2g
