/**
 * src/io.cpp
 * ==========
 * GFA v1.1 reader and GAF/VCF-inspired output writers.
 *
 * Phase 3 Task P3-8: GFA reader — COMPLETE.
 * Phase 3 Task P3-9: Output writer — COMPLETE.
 *
 * GFA line types handled:
 *   H  — header; ignored (version tag checked for info only).
 *   S  — segment: S <name> <sequence> [optional tags]
 *   L  — link:    L <u> <+/-> <v> <+/-> <overlap CIGAR>
 *   P  — path:    P <name> <seg1+,seg2+,...> <overlaps>
 *   W  — walk (GFA 1.1): W <sample> <hap_idx> <seq_id> <start> <end> <walk>
 *        Walk lines are parsed and stored as paths named "<sample>#<hap_idx>".
 *   #  — comment; ignored.
 *   Everything else: silently skipped.
 *
 * Orientation:
 *   GFA allows + (forward) and - (reverse-complement) orientations on both
 *   L and P lines.  Reverse-complement of a node sequence is computed on the
 *   fly when a '-' strand segment is encountered in a path.  The VariationGraph
 *   stores only forward-strand node sequences; the reader inserts a new
 *   reverse-complement node for each '-' occurrence (deduplicated by name
 *   "<name>_rc").  This is the same strategy used by vg and minigraph-cactus.
 *
 * Error handling:
 *   Malformed lines throw std::runtime_error with a message that includes the
 *   line number and the offending line content.
 */

#include <g2g/io.hpp>
#include <g2g/sequence.hpp>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace g2g {

// ── Internal helpers ──────────────────────────────────────────

/** Reverse-complement of a DNA sequence. */
static std::string rev_comp(std::string_view seq) {
    auto complement = [](char c) -> char {
        switch (c) {
            case 'A': case 'a': return 'T';
            case 'T': case 't': return 'A';
            case 'C': case 'c': return 'G';
            case 'G': case 'g': return 'C';
            default:            return 'N';
        }
    };
    std::string rc(seq.rbegin(), seq.rend());
    for (char& ch : rc) ch = complement(ch);
    return rc;
}

/** Split a string_view on a single delimiter character. */
static std::vector<std::string_view> split(std::string_view s, char delim) {
    std::vector<std::string_view> parts;
    size_t start = 0;
    for (size_t i = 0; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == delim) {
            if (i > start)
                parts.push_back(s.substr(start, i - start));
            start = i + 1;
        }
    }
    return parts;
}

/** Strip a trailing '\r' (for Windows-style line endings). */
static std::string_view strip_cr(std::string_view s) {
    if (!s.empty() && s.back() == '\r') return s.substr(0, s.size() - 1);
    return s;
}

// ── GFA parser ────────────────────────────────────────────────

/**
 * Core parsing logic shared by both read_gfa() overloads.
 *
 * Two-pass approach:
 *   Pass 1: parse S lines → register all nodes.
 *   Pass 2: parse L and P/W lines → add edges and paths.
 *
 * Why two passes?  GFA does not require S lines to appear before the P/L
 * lines that reference them.  Minimap2 and minigraph-cactus both produce
 * files where this holds, but the spec does not guarantee it.
 */
static void parse_gfa(std::istream& in,
                      const std::string& graph_name,
                      VariationGraph& vg) {
    vg.name = graph_name;

    // Buffer all lines so we can do two passes without rewinding
    // (istream may not be seekable, e.g. a pipe).
    std::vector<std::string> lines;
    {
        std::string line;
        while (std::getline(in, line)) lines.push_back(std::move(line));
    }

    // ── Pass 1: S lines → nodes ───────────────────────────────
    int lineno = 0;
    for (const auto& raw_line : lines) {
        ++lineno;
        std::string_view line = strip_cr(raw_line);
        if (line.empty() || line[0] == '#' || line[0] == 'H') continue;
        if (line[0] != 'S') continue;

        // S <name> <sequence> [optional tags...]
        auto parts = split(line, '\t');
        if (parts.size() < 3) {
            throw std::runtime_error(
                "GFA parse error at line " + std::to_string(lineno) +
                ": S line requires at least 3 fields: " + std::string(line));
        }
        std::string name(parts[1]);
        std::string_view seq = parts[2];

        // GFA allows '*' as a placeholder for an absent sequence.
        // We store it as empty; the aligner handles empty sequences.
        if (seq == "*") seq = "";

        vg.add_node(name, seq);
    }

    // ── Pass 2: L and P/W lines ───────────────────────────────
    lineno = 0;
    for (const auto& raw_line : lines) {
        ++lineno;
        std::string_view line = strip_cr(raw_line);
        if (line.empty() || line[0] == '#' || line[0] == 'H' || line[0] == 'S')
            continue;

        char type = line[0];

        // ── L line: L <u_name> <u_orient> <v_name> <v_orient> <overlap> ──
        if (type == 'L') {
            auto parts = split(line, '\t');
            if (parts.size() < 5) {
                throw std::runtime_error(
                    "GFA parse error at line " + std::to_string(lineno) +
                    ": L line requires 5 fields: " + std::string(line));
            }
            std::string u_name(parts[1]);
            char u_orient = parts[2].empty() ? '+' : parts[2][0];
            std::string v_name(parts[3]);
            char v_orient = parts[4].empty() ? '+' : parts[4][0];

            // Resolve orientations: for '-' strand, use/create the rc node.
            auto resolve = [&](const std::string& name, char orient) -> NodeId {
                if (orient == '+') {
                    auto it = vg.id_map.find(name);
                    if (it == vg.id_map.end())
                        throw std::runtime_error(
                            "GFA parse error at line " + std::to_string(lineno) +
                            ": unknown node '" + name + "'");
                    return it->second;
                } else {
                    // '-' orientation: use/create a reverse-complement node
                    std::string rc_name = name + "_rc";
                    auto it = vg.id_map.find(rc_name);
                    if (it != vg.id_map.end()) return it->second;
                    // Create rc node
                    auto src_it = vg.id_map.find(name);
                    if (src_it == vg.id_map.end())
                        throw std::runtime_error(
                            "GFA parse error at line " + std::to_string(lineno) +
                            ": unknown node '" + name + "'");
                    std::string rc_seq = rev_comp(vg.seq(src_it->second));
                    return vg.add_node(rc_name, rc_seq);
                }
            };

            NodeId u = resolve(u_name, u_orient);
            NodeId v = resolve(v_name, v_orient);
            vg.add_edge(u, v);
            continue;
        }

        // ── P line: P <path_name> <seg1+,seg2-,...> <overlaps> ───────────
        if (type == 'P') {
            auto parts = split(line, '\t');
            if (parts.size() < 3) {
                throw std::runtime_error(
                    "GFA parse error at line " + std::to_string(lineno) +
                    ": P line requires 3 fields: " + std::string(line));
            }
            std::string path_name(parts[1]);
            auto seg_parts = split(parts[2], ',');

            std::vector<NodeId> nodes;
            nodes.reserve(seg_parts.size());
            for (auto seg : seg_parts) {
                if (seg.empty()) continue;
                char orient = seg.back();
                if (orient == '+' || orient == '-')
                    seg = seg.substr(0, seg.size() - 1);
                else
                    orient = '+';

                std::string seg_name(seg);
                if (orient == '-') seg_name += "_rc";
                auto it = vg.id_map.find(seg_name);
                if (it == vg.id_map.end()) {
                    // rc node may not exist yet if there was no L line for it
                    if (orient == '-') {
                        auto src_it = vg.id_map.find(std::string(seg));
                        if (src_it == vg.id_map.end())
                            throw std::runtime_error(
                                "GFA parse error at line " + std::to_string(lineno) +
                                ": unknown segment '" + seg_name + "'");
                        std::string rc_seq = rev_comp(vg.seq(src_it->second));
                        NodeId rc_id = vg.add_node(seg_name, rc_seq);
                        nodes.push_back(rc_id);
                        continue;
                    }
                    throw std::runtime_error(
                        "GFA parse error at line " + std::to_string(lineno) +
                        ": unknown segment '" + seg_name + "'");
                }
                nodes.push_back(it->second);
            }
            if (!nodes.empty()) vg.add_path(path_name, nodes);
            continue;
        }

        // ── W line (GFA 1.1 walk): W <sample> <hap_idx> <seq_id> <start> <end> <walk> ──
        if (type == 'W') {
            auto parts = split(line, '\t');
            if (parts.size() < 7) {
                throw std::runtime_error(
                    "GFA parse error at line " + std::to_string(lineno) +
                    ": W line requires 7 fields: " + std::string(line));
            }
            // Path name: "<sample>#<hap_idx>#<seq_id>:<start>-<end>"
            std::string path_name = std::string(parts[1]) + "#" +
                                    std::string(parts[2]) + "#" +
                                    std::string(parts[3]) + ":" +
                                    std::string(parts[4]) + "-" +
                                    std::string(parts[5]);
            std::string_view walk = parts[6];

            // Walk format: >seg1<seg2>seg3  (> = forward, < = reverse)
            std::vector<NodeId> nodes;
            size_t i = 0;
            while (i < walk.size()) {
                char orient_ch = walk[i];
                if (orient_ch != '>' && orient_ch != '<') {
                    // Tolerate walks that start without orientation marker
                    ++i; continue;
                }
                char orient = (orient_ch == '>') ? '+' : '-';
                ++i;
                size_t j = i;
                while (j < walk.size() && walk[j] != '>' && walk[j] != '<') ++j;
                std::string seg_name(walk.substr(i, j - i));
                if (orient == '-') seg_name += "_rc";
                auto it = vg.id_map.find(seg_name);
                if (it == vg.id_map.end()) {
                    if (orient == '-') {
                        std::string base_name(walk.substr(i, j - i));
                        auto src = vg.id_map.find(base_name);
                        if (src == vg.id_map.end())
                            throw std::runtime_error(
                                "GFA parse error at line " + std::to_string(lineno) +
                                ": unknown segment '" + base_name + "' in W line");
                        std::string rc_seq = rev_comp(vg.seq(src->second));
                        NodeId rc_id = vg.add_node(seg_name, rc_seq);
                        nodes.push_back(rc_id);
                    } else {
                        throw std::runtime_error(
                            "GFA parse error at line " + std::to_string(lineno) +
                            ": unknown segment '" + seg_name + "' in W line");
                    }
                } else {
                    nodes.push_back(it->second);
                }
                i = j;
            }
            if (!nodes.empty()) vg.add_path(path_name, nodes);
            continue;
        }

        // All other line types (C, J, etc.): silently skip.
    }

    vg.finalise();
}

// ── Public API ────────────────────────────────────────────────

VariationGraph read_gfa(const std::filesystem::path& path) {
    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open GFA file: " + path.string());
    VariationGraph vg;
    parse_gfa(f, path.stem().string(), vg);
    return vg;
}

VariationGraph read_gfa(std::istream& in, const std::string& name) {
    VariationGraph vg;
    parse_gfa(in, name, vg);
    return vg;
}

// ── Output writers ────────────────────────────────────────────

void write_alignment_gaf(std::ostream& out,
                          const AlignmentResult& result,
                          const VariationGraph& ga,
                          const VariationGraph& gb,
                          const SnarlTree& ta,
                          const SnarlTree& tb) {
    // Header
    out << "# G2G alignment GAF\n";
    out << "# graph_a\t" << ga.name << "\n";
    out << "# graph_b\t" << gb.name << "\n";
    out << "# total_score\t" << result.total_score << "\n";
    out << "# anchors\t" << result.anchors.size() << "\n";
    out << "# snarl_alignments\t" << result.snarl_alignments.size() << "\n";

    // Anchor lines
    for (const auto& anc : result.anchors) {
        std::string_view na = (anc.node_a < ga.node_names.size())
                              ? std::string_view(ga.node_names[anc.node_a]) : "?";
        std::string_view nb = (anc.node_b < gb.node_names.size())
                              ? std::string_view(gb.node_names[anc.node_b]) : "?";
        out << "anchor\t" << na << "\t" << nb << "\t" << anc.sim << "\n";
    }

    // Snarl alignment lines
    for (const auto& sa : result.snarl_alignments) {
        // Source/sink names for context
        auto snarl_name = [](const VariationGraph& vg, const SnarlTree& tree,
                              uint32_t idx) -> std::string {
            if (idx >= tree.snarls.size()) return "?";
            const Snarl& s = tree.snarls[idx];
            std::string src = (s.source < vg.node_names.size())
                              ? vg.node_names[s.source] : std::to_string(s.source);
            std::string snk = (s.sink   < vg.node_names.size())
                              ? vg.node_names[s.sink]   : std::to_string(s.sink);
            return src + ">" + snk;
        };

        out << "snarl"
            << "\t" << snarl_name(ga, ta, sa.snarl_a_idx)
            << "\t" << snarl_name(gb, tb, sa.snarl_b_idx)
            << "\t" << sa.score
            << "\t" << sa.allele_matches.size();

        for (const auto& am : sa.allele_matches)
            out << "\t" << am.idx_a << ":" << am.idx_b << ":" << am.sim;
        out << "\n";
    }
}

void write_delta_vcf(std::ostream& out,
                      const AlignmentResult& result) {
    out << "##fileformat=G2GdeltaVCF-0.2\n";
    out << "##source=G2G graph-to-graph diff\n";
    out << "##note=One record per allele present in only one of the two graphs. "
           "With --ref: CHROM is the reference path, POS is the 0-BASED reference "
           "offset of the site (standard VCF is 1-based: add 1 to compare), REF/ALT "
           "are the atomic reference/alternate sequences (REF '.' = insertion, "
           "ALT '.' = deletion). Without --ref: CHROM is the snarl id, POS is 0 "
           "(graph space, no reference coordinate), REF is '.' and ALT is the "
           "sequence unique to one graph.\n";
    out << "##INFO=<ID=GRAPH,Number=1,Type=String,Description=\"Graph the allele "
           "is unique to (A=first input, B=second input)\">\n";
    out << "##INFO=<ID=SNARL,Number=1,Type=String,Description=\"Snarl region, "
           "source>sink\">\n";
    out << "##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"ALT sequence "
           "length\">\n";
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    // Emit every allele present in only one of the two graphs.  The sequences
    // come straight from AlignmentResult's delta lists (DeltaAllele::sequence),
    // which build_deltas() already populated — so the actual differing bases
    // are written out, not a placeholder.
    auto emit = [&](const std::vector<DeltaAllele>& deltas, const char* graph) {
        std::size_t i = 0;
        for (const auto& d : deltas) {
            // CHROM/POS: real reference coordinates when the diff supplied them
            // (reference-coordinate mode), else the snarl id with POS 0.
            if (d.pos >= 0) out << d.chrom << "\t" << d.pos;
            else            out << d.snarl_id << "\t0";
            out << "\t" << graph << '_' << i++         // ID
                << "\t" << (d.ref.empty() ? "." : d.ref)            // REF
                << "\t" << (d.sequence.empty() ? "." : d.sequence)  // ALT
                << "\t.\tPASS"
                << "\tGRAPH=" << graph
                << ";SNARL="  << d.snarl_id
                << ";LEN="    << d.sequence.size()
                << "\n";
        }
    };
    emit(result.a_only_alleles, "A");
    emit(result.b_only_alleles, "B");
}

} // namespace g2g
