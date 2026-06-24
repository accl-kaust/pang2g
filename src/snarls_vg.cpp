/**
 * src/snarls_vg.cpp
 * =================
 * Optional vg-backed snarl boundary provider.
 *
 * The built-in decomposer (decompose.cpp) finds snarl boundaries with a
 * dominator approximation that is correct only on DAGs and only reports
 * top-level snarls.  This module instead asks the `vg` toolkit's snarl finder
 * for the boundaries, which is the reference implementation used across the
 * pangenome ecosystem.
 *
 * Rather than link libvgio/libbdsg (a heavy Protobuf dependency), we shell out:
 *
 *     vg snarls -n <graph.gfa> | vg view -Rj -
 *
 *   -n   emit boundaries in named (GFA segment) space, so they map back to
 *        G2G node names directly.
 *   -Rj  read the snarl protobuf stream and print one JSON object per snarl.
 *
 * Each JSON line looks like:
 *   {"end": {"name": "m"}, "start": {"name": "a"}, "type": 1, ...}
 *
 * We parse the start/end segment names and translate them to G2G NodeIds.
 * Returns std::nullopt (so the caller falls back to the built-in decomposer)
 * if vg cannot be located or produces no usable output.
 *
 * Enabled at runtime with G2G_SNARLS=vg.  The vg binary is located via
 * G2G_VG (explicit path) or the PATH ("vg").
 */

#include <g2g/decompose.hpp>
#include <g2g/graph.hpp>

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include <unistd.h>   // getpid

namespace g2g {
namespace {

namespace fs = std::filesystem;

// Locate the vg binary: $G2G_VG if set, else "vg" (resolved via PATH).
std::string locate_vg() {
    if (const char* p = std::getenv("G2G_VG"); p && *p) return p;
    return "vg";
}

// Write the graph topology (H + S + L lines) to a GFA file.  Paths are not
// needed: vg's snarl finder works purely on topology.  Segment names are the
// G2G node names so the -n output maps straight back via id_map.
bool write_topology_gfa(const VariationGraph& vg, const fs::path& out) {
    std::ofstream f(out);
    if (!f) return false;
    f << "H\tVN:Z:1.1\n";
    for (NodeId n = 0; n < vg.num_nodes(); ++n) {
        std::string_view seq = vg.seq(n);
        const std::string& name = vg.node_names[n];
        f << "S\t" << name << "\t";
        f.write(seq.data(), static_cast<std::streamsize>(seq.size()));
        f << "\n";
    }
    for (NodeId u = 0; u < vg.num_nodes(); ++u)
        for (NodeId v : vg.successors(u))
            f << "L\t" << vg.node_names[u] << "\t+\t"
              << vg.node_names[v] << "\t+\t0M\n";
    return static_cast<bool>(f);
}

// Run a shell command, capturing stdout.  Returns nullopt if the pipeline
// could not be started.
std::optional<std::string> run_capture(const std::string& cmd) {
    std::FILE* pipe = ::popen(cmd.c_str(), "r");
    if (!pipe) return std::nullopt;
    std::string out;
    std::array<char, 8192> buf;
    std::size_t n;
    while ((n = std::fread(buf.data(), 1, buf.size(), pipe)) > 0)
        out.append(buf.data(), n);
    ::pclose(pipe);
    return out;
}

// Extract the segment name from the object that follows `key` (e.g. "start"):
//   ..."start": {"name": "a"}...        -> "a"
//   ..."end": {"backward": true, "name": "x"}... -> "x"
// Returns the first quoted value after the next "name" token past `key`.
std::optional<std::string> name_after(const std::string& s, const char* key) {
    std::size_t kp = s.find(key);
    if (kp == std::string::npos) return std::nullopt;
    std::size_t np = s.find("\"name\"", kp);
    if (np == std::string::npos) return std::nullopt;
    std::size_t colon = s.find(':', np);
    if (colon == std::string::npos) return std::nullopt;
    std::size_t q1 = s.find('"', colon);
    if (q1 == std::string::npos) return std::nullopt;
    std::size_t q2 = s.find('"', q1 + 1);
    if (q2 == std::string::npos) return std::nullopt;
    return s.substr(q1 + 1, q2 - q1 - 1);
}

}  // namespace

std::optional<std::vector<std::pair<NodeId, NodeId>>>
vg_snarl_boundaries(const VariationGraph& vg) {
    if (vg.num_nodes() == 0)
        return std::vector<std::pair<NodeId, NodeId>>{};

    // Unique temp GFA in the system temp dir.
    std::error_code ec;
    fs::path tmp = fs::temp_directory_path(ec) /
        ("g2g_vg_" + std::to_string(::getpid()) + "_" +
         std::to_string(reinterpret_cast<std::uintptr_t>(&vg)) + ".gfa");
    if (ec) return std::nullopt;
    if (!write_topology_gfa(vg, tmp)) { fs::remove(tmp, ec); return std::nullopt; }

    const std::string vgbin = locate_vg();
    // vg snarls -n <gfa> | vg view -Rj -      (stderr silenced; we detect
    // failure via empty/again output and fall back).
    std::string cmd = vgbin + " snarls -n " + "\"" + tmp.string() + "\""
                    + " 2>/dev/null | " + vgbin + " view -Rj - 2>/dev/null";
    std::optional<std::string> out = run_capture(cmd);
    fs::remove(tmp, ec);

    if (!out || out->empty()) return std::nullopt;

    std::vector<std::pair<NodeId, NodeId>> boundaries;
    std::size_t pos = 0;
    while (pos < out->size()) {
        std::size_t nl = out->find('\n', pos);
        std::string line = out->substr(pos, nl == std::string::npos
                                                ? std::string::npos : nl - pos);
        pos = (nl == std::string::npos) ? out->size() : nl + 1;
        if (line.find("\"start\"") == std::string::npos) continue;

        auto sn = name_after(line, "\"start\"");
        auto en = name_after(line, "\"end\"");
        if (!sn || !en) continue;

        auto its = vg.id_map.find(*sn);
        auto ite = vg.id_map.find(*en);
        if (its == vg.id_map.end() || ite == vg.id_map.end()) continue;
        boundaries.emplace_back(its->second, ite->second);
    }

    // No parsed boundaries despite output → treat as failure (fall back).
    if (boundaries.empty()) return std::nullopt;
    return boundaries;
}

}  // namespace g2g
