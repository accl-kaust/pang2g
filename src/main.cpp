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

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>

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
        "  --sequential       Use single-threaded aligner\n"
        "  --snarls <backend> Snarl decomposition: 'builtin' (default, DAG-only)\n"
        "                     or 'vg' (use the vg toolkit; needs vg on PATH or\n"
        "                     $G2G_VG; falls back to builtin if unavailable)\n"
        "  --verbose-timing   Print [timing] (io/index/decompose/dp) and, for\n"
        "                     the parallel runtime, [runstats] (tasks/steals/\n"
        "                     imbalance) lines to stderr\n"
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
        else if (flag == "--sequential") a.sequential = true;
        else if (flag == "--snarls") {
            a.snarls = need_next();
            if (a.snarls != "vg" && a.snarls != "builtin") {
                std::cerr << "ERROR: --snarls expects 'vg' or 'builtin'\n";
                std::exit(1);
            }
        }
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

    // ── Run alignment ─────────────────────────────────────────
    g2g::AlignmentResult result;
    g2g::RunStats run_stats{};   // .phases filled in both paths; scheduler
                                 // fields (steals/imbalance/...) only parallel
    bool have_runstats = false;
    try {
        if (args.sequential) {
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

    // ── Verbose timing / runtime stats ────────────────────────
    if (args.verbose_timing) {
        std::cerr << std::fixed << std::setprecision(6)
                  << "[timing] io="    << io_seconds
                  << " index="         << run_stats.phases.index_seconds
                  << " decompose="     << run_stats.phases.decompose_seconds
                  << " dp="            << run_stats.phases.dp_seconds
                  << "\n";
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

    // ── Decompose for writers (need SnarlTree with node names) ─
    g2g::SnarlTree ta = g2g::SnarlDecomposer(ga).decompose();
    g2g::SnarlTree tb = g2g::SnarlDecomposer(gb).decompose();

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
