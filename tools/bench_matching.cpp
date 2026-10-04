/**
 * tools/bench_matching.cpp
 * ========================
 * Measure the per-snarl-pair bipartite matching cost as a function of the
 * allele count k, for the exact LAPJV solver and for brute-force enumeration
 * (the pre-LAPJV baseline, feasible only for small k).  Produces the CSV that
 * make_figures.py turns into fig_matching (columns: k,brute_ms,lapjv_ms).
 *
 * Usage:
 *   bench_matching [--k-max 50] [--brute-k-max 8] [--seq-len 12] [--reps 20]
 *                  [--seed 42] > results/matching_cost_vs_k.csv
 *
 * Each row times `reps` independent snarl pairs with k random alleles per
 * side (half of the B alleles are copies of A alleles, so the matching is
 * non-trivial but has a clear optimum), and reports the mean wall time per
 * pair in milliseconds.  Brute force is skipped (blank) above --brute-k-max.
 */

#include <g2g/matching.hpp>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <string_view>
#include <vector>

using namespace g2g;

static std::string random_dna(int len, std::mt19937& rng) {
    static const char B[] = "ACGT";
    std::uniform_int_distribution<int> d(0, 3);
    std::string s(static_cast<std::size_t>(len), 'A');
    for (char& c : s) c = B[d(rng)];
    return s;
}

template <class F>
static double time_ms(F&& f, int reps) {
    using Clock = std::chrono::steady_clock;
    auto t0 = Clock::now();
    for (int r = 0; r < reps; ++r) f(r);
    auto t1 = Clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count() / reps;
}

int main(int argc, char** argv) {
    int k_max = 50, brute_k_max = 8, seq_len = 12, reps = 20;
    unsigned seed = 42;
    for (int i = 1; i < argc; ++i) {
        auto next = [&](int& out) { if (i + 1 < argc) out = std::atoi(argv[++i]); };
        if (!std::strcmp(argv[i], "--k-max")) next(k_max);
        else if (!std::strcmp(argv[i], "--brute-k-max")) next(brute_k_max);
        else if (!std::strcmp(argv[i], "--seq-len")) next(seq_len);
        else if (!std::strcmp(argv[i], "--reps")) next(reps);
        else if (!std::strcmp(argv[i], "--seed")) { int s = 42; next(s); seed = static_cast<unsigned>(s); }
    }

    std::printf("k,brute_ms,lapjv_ms\n");
    for (int k = 2; k <= k_max; ++k) {
        std::mt19937 rng(seed + static_cast<unsigned>(k));
        // Pre-generate `reps` instances so timing excludes generation.
        std::vector<std::vector<std::string>> A(reps), Bv(reps);
        for (int r = 0; r < reps; ++r) {
            for (int i = 0; i < k; ++i) A[r].push_back(random_dna(seq_len, rng));
            for (int i = 0; i < k; ++i)
                Bv[r].push_back(i % 2 == 0 ? A[r][static_cast<std::size_t>(i)]
                                           : random_dna(seq_len, rng));
        }
        auto views = [](const std::vector<std::string>& v) {
            std::vector<std::string_view> out;
            out.reserve(v.size());
            for (const auto& s : v) out.emplace_back(s);
            return out;
        };
        std::vector<std::vector<std::string_view>> Av(reps), Bvv(reps);
        for (int r = 0; r < reps; ++r) { Av[r] = views(A[r]); Bvv[r] = views(Bv[r]); }

        MatchParams mp;
        mp.lapjv_limit = 1 << 30;   // force LAPJV at every k (no greedy fallback)

        volatile float sink = 0.f;
        double lapjv_ms = time_ms([&](int r) {
            MatchResult m = hungarian_match(Av[r], Bvv[r], mp);
            sink = sink + m.score;
        }, reps);

        if (k <= brute_k_max) {
            int breps = k >= 8 ? 1 : (k >= 7 ? 2 : reps);
            double brute_ms = time_ms([&](int r) {
                MatchResult m = brute_match(Av[r % reps], Bvv[r % reps], mp);
                sink = sink + m.score;
            }, breps);
            std::printf("%d,%.6f,%.6f\n", k, brute_ms, lapjv_ms);
        } else {
            std::printf("%d,,%.6f\n", k, lapjv_ms);
        }
        std::fflush(stdout);
    }
    return 0;
}
