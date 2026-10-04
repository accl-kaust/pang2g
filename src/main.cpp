/**
 * src/main.cpp
 * ============
 * G2G command-line interface.
 *
 * Usage:
 *   g2g <graph_a.gfa> <graph_b.gfa> [options]
 *
 * Options:
 *   --gaf <file>       Write alignment output in GAF-like format (default: stdout)
 *   --vcf <file>       Write A-only delta alleles in VCF-inspired format
 *   --threads <n>      Number of threads (default: all hardware threads)
 *   --alpha <f>        W_seq weight       (default: 1.0)
 *   --beta  <f>        W_topo weight      (default: 0.5)
 *   --gamma <f>        W_struct weight    (default: 2.0)
 *   --gap   <f>        Gap penalty        (default: 0.5)
 *   --sequential       Use single-threaded aligner (skip Phase 4 runtime)
 *   -h / --help        Print this help
 *
 * Exit codes:
 *   0  success
 *   1  bad arguments or I/O error
 *   2  alignment produced no snarl pairs (may indicate malformed input)
 */

#include <g2g/align.hpp>
#include <g2g/decompose.hpp>
#include <g2g/io.hpp>
#include <g2g/parallel.hpp>
#include <g2g/refdiff.hpp>
#include <g2g/sequence.hpp>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>

#include <sys/resource.h>   // getrusage — self-reported peak RSS (no GNU time dep)

namespace fs = std::filesystem;

// ── Argument parsing ──────────────────────────────────────────

struct Args {
    fs::path gfa_a;
    fs::path gfa_b;
    fs::path gaf_out;       // empty = stdout
    fs::path vcf_out;       // empty = not written
    int      threads    = 0;  // 0 = auto
    bool     sequential = false;
    bool     verbose_timing = false;
    std::string snarls;     // "vg" | "builtin" | "" (default = builtin)
    std::string ref;        // reference path name → reference-coordinate diff
    fs::path dump_snarls;   // write the decomposed snarls of both graphs (debug)
    g2g::AlignerParams params;
};

static void print_help(const char* prog) {
    std::cout <<
        "Usage: " << prog << " <graph_a.gfa> <graph_b.gfa> [options]\n"
        "\n"
        "Options:\n"
        "  --gaf <file>       GAF-like alignment output (default: stdout)\n"
        "  --vcf <file>       VCF-inspired delta output (A-only alleles)\n"
        "  --threads <n>      Worker threads (default: hardware_concurrency)\n"
        "  --alpha <f>        W_seq weight   [default: 1.0]\n"
        "  --beta  <f>        W_topo weight  [default: 0.5]\n"
        "  --gamma <f>        W_struct weight[default: 2.0]\n"
        "  --gap   <f>        Gap penalty    [default: 0.5]\n"
        "  --ref <path>       Reference-coordinate diff: match snarls by their\n"
        "                     position along the embedded reference path <path>\n"
        "                     (robust to differing node segmentation; emits real\n"
        "                     coordinates). Both graphs must contain that path.\n"
        "  --sequential       Use single-threaded aligner (ignored with --ref)\n"
        "  --pairing <mode>   Snarl pairing for the graph-to-graph aligner:\n"
        "                     'boundary' (default): identical source/sink boundary\n"
        "                     sequences, repeats consumed positionally, O(n);\n"
        "                     'interval': anchor-guided, order-preserving DP between\n"
        "                     consecutive shared anchors (contains pairing errors\n"
        "                     to one anchor interval). Ignored with --ref.\n"
        "  --snarls <backend> Snarl decomposition: 'builtin' (default, DAG-only)\n"
        "                     or 'vg' (use the vg toolkit; needs vg on PATH or\n"
        "                     $G2G_VG; falls back to builtin if unavailable)\n"
        "  --dump-snarls <f>  Write every decomposed snarl of both graphs to <f>\n"
        "                     (graph, source>sink, k, allele sequences) for debugging\n"
        "  --verbose-timing   Print [timing] (io/decompose/index/dp/wall) and\n"
        "                     [mem] (peak_rss_mb) lines, plus [runstats] (tasks/\n"
        "                     steals/imbalance) for the parallel runtime, to stderr\n"
        "  -h, --help         Show this help\n";
}

static Args parse_args(int argc, char** argv) {
    if (argc < 3) {
        print_help(argv[0]);
        std::exit(1);
    }

    Args a;
    a.gfa_a = argv[1];
    a.gfa_b = argv[2];

    for (int i = 3; i < argc; ++i) {
        std::string flag = argv[i];
        auto need_next = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "ERROR: " << flag << " requires an argument\n";
                std::exit(1);
            }
            return argv[++i];
        };

        if (flag == "--gaf")        a.gaf_out   = need_next();
        else if (flag == "--vcf")   a.vcf_out   = need_next();
        else if (flag == "--threads") a.threads  = std::stoi(need_next());
        else if (flag == "--alpha") a.params.alpha = std::stof(need_next());
        else if (flag == "--beta")  a.params.beta  = std::stof(need_next());
        else if (flag == "--gamma") a.params.gamma = std::stof(need_next());
        else if (flag == "--gap")   a.params.gap   = std::stof(need_next());
        else if (flag == "--ref") a.ref = need_next();
        else if (flag == "--sequential") a.sequential = true;
        else if (flag == "--pairing") {
            std::string m = need_next();
            if (m == "boundary")      a.params.pairing = g2g::PairingMode::Boundary;
            else if (m == "interval") a.params.pairing = g2g::PairingMode::Interval;
            else { std::cerr << "ERROR: --pairing expects 'boundary' or 'interval'\n"; std::exit(1); }
        }
        else if (flag == "--snarls") {
            a.snarls = need_next();
            if (a.snarls != "vg" && a.snarls != "builtin") {
                std::cerr << "ERROR: --snarls expects 'vg' or 'builtin'\n";
                std::exit(1);
            }
        }
        else if (flag == "--dump-snarls") a.dump_snarls = need_next();
        else if (flag == "--verbose-timing") a.verbose_timing = true;
        else if (flag == "-h" || flag == "--help") { print_help(argv[0]); std::exit(0); }
        else { std::cerr << "Unknown option: " << flag << "\n"; std::exit(1); }
    }

    for (const auto& p : { a.gfa_a, a.gfa_b }) {
        if (!fs::exists(p)) {
            std::cerr << "ERROR: file not found: " << p << "\n";
            std::exit(1);
        }
    }
    return a;
}

// ── Main ──────────────────────────────────────────────────────

int main(int argc, char** argv) {
    Args args = parse_args(argc, argv);

    // --snarls drives the decomposer back-end via the environment (read inside
    // SnarlDecomposer::decompose()); the explicit flag overrides any inherited
    // G2G_SNARLS.
    if (!args.snarls.empty())
        ::setenv("G2G_SNARLS", args.snarls.c_str(), /*overwrite=*/1);

    using Clock = std::chrono::steady_clock;
    const auto wall_t0 = Clock::now();

    // ── Load graphs ───────────────────────────────────────────
    g2g::VariationGraph ga, gb;
    double io_seconds = 0.0;
    try {
        auto io_t0 = Clock::now();
        std::cerr << "[g2g] reading " << args.gfa_a << " ...\n";
        ga = g2g::read_gfa(args.gfa_a);
        std::cerr << "[g2g] reading " << args.gfa_b << " ...\n";
        gb = g2g::read_gfa(args.gfa_b);
        io_seconds = std::chrono::duration<double>(Clock::now() - io_t0).count();
    } catch (const std::exception& e) {
        std::cerr << "ERROR reading GFA: " << e.what() << "\n";
        return 1;
    }
    std::cerr << "[g2g] GA: " << ga.num_nodes() << " nodes"
              << "  GB: " << gb.num_nodes() << " nodes\n";

    // ── Decompose both graphs (shared by the diff and the writers) ─
    g2g::SnarlTree ta, tb;
    double decompose_seconds = 0.0;
    try {
        auto dec_t0 = Clock::now();
        ta = g2g::SnarlDecomposer(ga).decompose();
        tb = g2g::SnarlDecomposer(gb).decompose();
        decompose_seconds = std::chrono::duration<double>(Clock::now() - dec_t0).count();
    } catch (const std::exception& e) {
        std::cerr << "ERROR during snarl decomposition: " << e.what() << "\n";
        return 1;
    }

    // ── Optional snarl dump (debugging aid) ───────────────────
    if (!args.dump_snarls.empty()) {
        std::ofstream f(args.dump_snarls);
        if (!f) { std::cerr << "ERROR: cannot open " << args.dump_snarls << "\n"; return 1; }
        f << "#graph\tsource>sink\tk\talleles\n";
        auto dump = [&](const char* tag, const g2g::VariationGraph& g, const g2g::SnarlTree& t) {
            std::vector<std::string> seqs;
            for (const auto& s : t.snarls) {
                g2g::allele_sequences(g, s, t, seqs);
                f << tag << "\t" << g.node_names[s.source] << ">" << g.node_names[s.sink]
                  << "\t" << s.allele_count() << "\t";
                for (std::size_t i = 0; i < seqs.size(); ++i)
                    f << (i ? "," : "") << (seqs[i].empty() ? "-" : seqs[i]);
                f << "\n";
            }
        };
        dump("A", ga, ta);
        dump("B", gb, tb);
        std::cerr << "[g2g] snarls written to " << args.dump_snarls << "\n";
    }

    // ── Run alignment / diff ──────────────────────────────────
    g2g::AlignmentResult result;
    g2g::RunStats run_stats{};   // .phases filled in both paths; scheduler
                                 // fields (steals/imbalance/...) only parallel
    bool have_runstats = false;
    try {
        if (!args.ref.empty()) {
            std::cerr << "[g2g] reference-coordinate diff along path '"
                      << args.ref << "'\n";
            auto rd = g2g::reference_diff(ga, gb, ta, tb, args.ref);
            if (!rd) {
                std::cerr << "ERROR: reference path '" << args.ref
                          << "' not found in both graphs\n";
                return 1;
            }
            result = std::move(*rd);
        } else if (args.sequential) {
            std::cerr << "[g2g] running sequential aligner\n";
            g2g::GraphAligner aligner(args.params);
            result = aligner.align(ga, gb);
            run_stats.phases = aligner.last_phase_times();
        } else {
            int nthreads = args.threads > 0
                ? args.threads
                : static_cast<int>(std::thread::hardware_concurrency());
            std::cerr << "[g2g] running parallel aligner (" << nthreads << " threads)\n";
            g2g::ParallelAligner aligner(args.params, nthreads);
            result = aligner.align(ga, gb);
            run_stats     = aligner.last_stats();
            have_runstats = true;
        }
    } catch (const std::exception& e) {
        std::cerr << "ERROR during alignment: " << e.what() << "\n";
        return 1;
    }

    std::cerr << "[g2g] total_score=" << result.total_score
              << "  snarl_pairs="     << result.snarl_alignments.size()
              << "  a_only="          << result.a_only_alleles.size()
              << "  b_only="          << result.b_only_alleles.size()
              << "\n";
    if (args.ref.empty()) {
        std::cerr << "[g2g] anchors=" << result.anchors.size()
                  << "  unpaired_snarls_a=" << result.unpaired_snarls_a
                  << "  unpaired_snarls_b=" << result.unpaired_snarls_b
                  << "  pairing=" << (args.params.pairing == g2g::PairingMode::Interval
                                        ? "interval" : "boundary")
                  << "\n";
    }

    // ── Verbose timing / runtime stats ────────────────────────
    if (args.verbose_timing) {
        const double wall_seconds =
            std::chrono::duration<double>(Clock::now() - wall_t0).count();
        // Self-reported peak resident memory — avoids depending on GNU
        // /usr/bin/time being installed on the host. On Linux ru_maxrss is KiB.
        struct rusage ru{};
        getrusage(RUSAGE_SELF, &ru);
        const double peak_rss_mb = static_cast<double>(ru.ru_maxrss) / 1024.0;

        std::cerr << std::fixed << std::setprecision(6)
                  << "[timing] io="    << io_seconds
                  << " decompose="     << decompose_seconds
                  << " index="         << run_stats.phases.index_seconds
                  << " dp="            << run_stats.phases.dp_seconds
                  << " wall="          << wall_seconds
                  << "\n";
        std::cerr << std::fixed << std::setprecision(1)
                  << "[mem] peak_rss_mb=" << peak_rss_mb << "\n";
        if (have_runstats) {
            std::cerr << std::setprecision(4)
                      << "[runstats] tasks=" << run_stats.tasks_executed
                      << " steals="          << run_stats.steals
                      << " imbalance="       << run_stats.imbalance
                      << " ready_peak="      << run_stats.tasks_ready_peak
                      << " threads="         << run_stats.num_threads
                      << "\n";
        }
        std::cerr.unsetf(std::ios::floatfield);
    }

    // (ta, tb already decomposed above and reused by the writers.)

    // ── Write GAF ─────────────────────────────────────────────
    if (args.gaf_out.empty()) {
        g2g::write_alignment_gaf(std::cout, result, ga, gb, ta, tb);
    } else {
        std::ofstream f(args.gaf_out);
        if (!f) { std::cerr << "ERROR: cannot open " << args.gaf_out << "\n"; return 1; }
        g2g::write_alignment_gaf(f, result, ga, gb, ta, tb);
        std::cerr << "[g2g] alignment written to " << args.gaf_out << "\n";
    }

    // ── Write VCF ─────────────────────────────────────────────
    if (!args.vcf_out.empty()) {
        std::ofstream f(args.vcf_out);
        if (!f) { std::cerr << "ERROR: cannot open " << args.vcf_out << "\n"; return 1; }
        g2g::write_delta_vcf(f, result);
        std::cerr << "[g2g] delta VCF written to " << args.vcf_out << "\n";
    }

    return result.snarl_alignments.empty() ? 2 : 0;
}
