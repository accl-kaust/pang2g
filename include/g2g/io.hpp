#pragma once
/**
 * include/g2g/io.hpp
 * ==================
 * GFA v1.1 reader/writer and GAF-like output format.
 *
 * GFA reader:
 *   - Parses S (segment), L (link), P (path) and W (walk) lines.
 *   - Stores sequences in flat VariationGraph::seq_data.
 *   - W (walk) lines are stored as paths named PanSN-style
 *     "<sample>#<hap>#<seqid>:<start>-<end>"; resolve them with --ref by
 *     sample or contig (see refdiff).
 *   - Handles both GFA 1.0 (P-lines) and GFA 1.1 (W-lines).
 *
 * Output format:
 *   write_alignment_gaf()  — one line per matched snarl pair:
 *       <snarl_a_id> <snarl_b_id> <score> <allele_matches...>
 *   write_delta_vcf()      — VCF-inspired file; one record per
 *       A-only allele (delta graph ΔA→B).
 *
 * Phase 3 task estimate: 2 weeks (§3.3, item 1 + item 6).
 */

#include <g2g/types.hpp>
#include <g2g/graph.hpp>
#include <g2g/decompose.hpp>   // add this line
#include <filesystem>
#include <iosfwd>
#include <string>

namespace g2g {

/**
 * Parse a GFA v1.1 file into a VariationGraph.
 * Throws std::runtime_error on malformed input.
 * Call vg.finalise() automatically before returning.
 */
VariationGraph read_gfa(const std::filesystem::path& path);

/**
 * Parse from an already-open stream (useful for gzip piping:
 *   process_gfa("zcat graph.gfa.gz | g2g ...") via a popen stream).
 */
VariationGraph read_gfa(std::istream& in, const std::string& name = "");

/**
 * Write an AlignmentResult to a GAF-like plain-text file.
 * Format (tab-separated):
 *   snarl_a   snarl_b   score   n_matches   [i:j:sim ...]
 */
void write_alignment_gaf(std::ostream& out,
                          const AlignmentResult& result,
                          const VariationGraph& ga,
                          const VariationGraph& gb,
                          const SnarlTree& ta,
                          const SnarlTree& tb);

/**
 * Write the delta graph in VCF-inspired format: one record per allele present
 * in only one of the two graphs.  ALT carries the actual unique sequence;
 * GRAPH=A marks alleles unique to the first input, GRAPH=B the second.
 * CHROM = snarl id; POS = 0 (graph space, no reference coordinate); REF = ".".
 */
void write_delta_vcf(std::ostream& out,
                      const AlignmentResult& result);

} // namespace g2g
