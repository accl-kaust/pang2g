/**
 * tests/test_pairing.cpp
 * ======================
 * Snarl pairing and delta construction.
 *
 *  * Unpaired snarls: every allele of a snarl that finds no partner must reach
 *    the delta (Definition 2), in both the sequential and parallel aligners.
 *  * Boundary-key cascade: with short, repeated boundary sequences, one snarl
 *    private to A shifts the positional (FIFO) pairing of every later snarl.
 *  * Anchor-interval pairing (PairingMode::Interval) contains that error to
 *    one interval and recovers the planted deltas exactly; with no anchors it
 *    degrades to one order-preserving DP over the whole graph, which is still
 *    exact on these small inputs.
 */

#include <gtest/gtest.h>
#include <g2g/align.hpp>
#include <g2g/parallel.hpp>
#include <g2g/graph.hpp>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

using namespace g2g;

namespace {

// Deterministic 12-mer for (slot, allele) — identical in both graphs.
std::string allele_seq(int slot, int a) {
    static const char B[] = "ACGT";
    std::string s(12, 'A');
    unsigned x = static_cast<unsigned>(slot * 131 + a * 17 + 7);
    for (char& c : s) { x = x * 1103515245u + 12345u; c = B[(x >> 16) & 3]; }
    return s;
}
std::string linker_seq(int i) {
    static const char B[] = "ACGT";
    std::string s(24, 'A');
    unsigned x = static_cast<unsigned>(i * 977 + 3);
    for (char& c : s) { x = x * 1664525u + 1013904223u; c = B[(x >> 16) & 3]; }
    return s;
}

/**
 * Slot kinds: 'S' shared (2 shared alleles + 1 private allele per graph),
 * 'A' A-private bubble (3 alleles, B has a plain edge), 'B' symmetric.
 * `boundary` is the sequence of EVERY boundary node (short + identical →
 * colliding keys; pass distinct sequences via `distinct=true` instead).
 * `linker_every` inserts a shared 24-bp linker node (anchor) after every
 * m-th slot.  Returns graphs plus the planted A-only / B-only sequences.
 */
struct Pair {
    VariationGraph ga, gb;
    std::vector<std::string> truth_a, truth_b;
    int shared_slots{0};
};

Pair build(const std::string& kinds, bool distinct, int linker_every) {
    Pair p;
    p.ga.name = "A"; p.gb.name = "B";
    auto bseq = [&](int i) { return distinct ? "GATTACA" + std::to_string(i) + "C" : std::string("A"); };
    NodeId sa = p.ga.add_node("b0", bseq(0)), sb = p.gb.add_node("b0", bseq(0));
    std::vector<NodeId> path_a{sa}, path_b{sb};
    int bidx = 1;
    for (int s = 0; s < static_cast<int>(kinds.size()); ++s) {
        NodeId ta = p.ga.add_node("b" + std::to_string(bidx), bseq(bidx));
        NodeId tb = p.gb.add_node("b" + std::to_string(bidx), bseq(bidx));
        ++bidx;
        char k = kinds[static_cast<std::size_t>(s)];
        NodeId first_a = kInvalidNode, first_b = kInvalidNode;
        if (k == 'S') {
            ++p.shared_slots;
            for (int a = 0; a < 2; ++a) {
                std::string q = allele_seq(s, a);
                NodeId na = p.ga.add_node("s" + std::to_string(s) + "sh" + std::to_string(a) + "a", q);
                NodeId nb = p.gb.add_node("s" + std::to_string(s) + "sh" + std::to_string(a) + "b", q);
                p.ga.add_edge(sa, na); p.ga.add_edge(na, ta);
                p.gb.add_edge(sb, nb); p.gb.add_edge(nb, tb);
                if (a == 0) { first_a = na; first_b = nb; }
            }
            std::string qa = allele_seq(s, 10), qb = allele_seq(s, 11);
            p.truth_a.push_back(qa); p.truth_b.push_back(qb);
            NodeId na = p.ga.add_node("s" + std::to_string(s) + "ao", qa);
            NodeId nb = p.gb.add_node("s" + std::to_string(s) + "bo", qb);
            p.ga.add_edge(sa, na); p.ga.add_edge(na, ta);
            p.gb.add_edge(sb, nb); p.gb.add_edge(nb, tb);
        } else {
            VariationGraph& own = (k == 'A') ? p.ga : p.gb;
            VariationGraph& oth = (k == 'A') ? p.gb : p.ga;
            std::vector<std::string>& tr = (k == 'A') ? p.truth_a : p.truth_b;
            NodeId src = (k == 'A') ? sa : sb, snk = (k == 'A') ? ta : tb;
            for (int a = 0; a < 3; ++a) {
                std::string q = allele_seq(s, 20 + a);
                tr.push_back(q);
                NodeId n = own.add_node("s" + std::to_string(s) + "p" + std::to_string(a), q);
                own.add_edge(src, n); own.add_edge(n, snk);
                if (a == 0) { if (k == 'A') first_a = n; else first_b = n; }
            }
            oth.add_edge((k == 'A') ? sb : sa, (k == 'A') ? tb : ta);
        }
        if (first_a != kInvalidNode) path_a.push_back(first_a);
        path_a.push_back(ta);
        if (first_b != kInvalidNode) path_b.push_back(first_b);
        path_b.push_back(tb);
        sa = ta; sb = tb;
        if (linker_every > 0 && (s + 1) % linker_every == 0 && s + 1 < static_cast<int>(kinds.size())) {
            std::string ls = linker_seq(s);
            NodeId la = p.ga.add_node("lnk" + std::to_string(s), ls);
            NodeId lb = p.gb.add_node("lnk" + std::to_string(s), ls);
            NodeId na = p.ga.add_node("b" + std::to_string(bidx), bseq(bidx));
            NodeId nb = p.gb.add_node("b" + std::to_string(bidx), bseq(bidx));
            ++bidx;
            p.ga.add_edge(sa, la); p.ga.add_edge(la, na);
            p.gb.add_edge(sb, lb); p.gb.add_edge(lb, nb);
            path_a.push_back(la); path_a.push_back(na);
            path_b.push_back(lb); path_b.push_back(nb);
            sa = na; sb = nb;
        }
    }
    p.ga.add_path("hapA", path_a);
    p.gb.add_path("hapB", path_b);
    p.ga.finalise(); p.gb.finalise();
    return p;
}

std::multiset<std::string> seqs(const std::vector<DeltaAllele>& d) {
    std::multiset<std::string> out;
    for (const auto& x : d) out.insert(x.sequence);
    return out;
}
std::multiset<std::string> ms(const std::vector<std::string>& v) {
    return std::multiset<std::string>(v.begin(), v.end());
}

AlignerParams params(PairingMode m) { AlignerParams p; p.pairing = m; return p; }

}  // namespace

// ── Unpaired snarls reach the delta ──────────────────────────────────────────

TEST(Deltas, UnpairedSnarlAllelesAreReported) {
    Pair p = build("SASS", /*distinct=*/true, /*linker_every=*/0);
    GraphAligner seq(params(PairingMode::Boundary));
    AlignmentResult r = seq.align(p.ga, p.gb);
    EXPECT_EQ(r.snarl_alignments.size(), 3u);          // the three shared slots
    EXPECT_EQ(r.unpaired_snarls_a, 1u);                // the A-private bubble
    EXPECT_EQ(r.unpaired_snarls_b, 0u);
    EXPECT_EQ(seqs(r.a_only_alleles), ms(p.truth_a));  // includes its 3 alleles
    EXPECT_EQ(seqs(r.b_only_alleles), ms(p.truth_b));

    ParallelAligner par(params(PairingMode::Boundary), 3);
    AlignmentResult q = par.align(p.ga, p.gb);
    EXPECT_EQ(seqs(q.a_only_alleles), seqs(r.a_only_alleles));
    EXPECT_EQ(seqs(q.b_only_alleles), seqs(r.b_only_alleles));
    EXPECT_EQ(q.unpaired_snarls_a, r.unpaired_snarls_a);
}

// ── Boundary-key cascade ─────────────────────────────────────────────────────

TEST(Pairing, BoundaryKeyCollisionCascades) {
    // 1-bp identical boundaries everywhere; slot 1 exists only in A.
    Pair p = build("SASSSSSS", /*distinct=*/false, /*linker_every=*/0);
    GraphAligner seq(params(PairingMode::Boundary));
    AlignmentResult r = seq.align(p.ga, p.gb);
    // Every key is identical, so A-snarl i is paired with B-snarl i: slot 0 is
    // right, slot 1 (private) takes B's slot 2, and every later pair is shifted.
    EXPECT_EQ(r.snarl_alignments.size(), 7u);
    EXPECT_EQ(r.unpaired_snarls_a, 1u);
    EXPECT_EQ(r.unpaired_snarls_b, 0u);
    // Mispaired snarls match no alleles, so their alleles all land in the
    // delta: many more delta alleles than were planted.
    EXPECT_GT(r.a_only_alleles.size(), p.truth_a.size());
    EXPECT_GT(r.b_only_alleles.size(), p.truth_b.size());
    EXPECT_NE(seqs(r.a_only_alleles), ms(p.truth_a));
}

// ── Anchor-interval pairing contains the error ───────────────────────────────

TEST(Pairing, IntervalPairingRecoversPlantedDeltas) {
    Pair p = build("SASSSSSS", /*distinct=*/false, /*linker_every=*/3);
    GraphAligner seq(params(PairingMode::Interval));
    AlignmentResult r = seq.align(p.ga, p.gb);
    EXPECT_GE(r.anchors.size(), 2u);                   // the shared linkers
    EXPECT_EQ(r.snarl_alignments.size(), static_cast<std::size_t>(p.shared_slots));
    EXPECT_EQ(r.unpaired_snarls_a, 1u);
    EXPECT_EQ(r.unpaired_snarls_b, 0u);
    for (const auto& aln : r.snarl_alignments)
        EXPECT_EQ(aln.allele_matches.size(), 2u) << "paired snarls share exactly 2 alleles";
    EXPECT_EQ(seqs(r.a_only_alleles), ms(p.truth_a));
    EXPECT_EQ(seqs(r.b_only_alleles), ms(p.truth_b));

    // Parallel runtime pairs identically.
    ParallelAligner par(params(PairingMode::Interval), 4);
    AlignmentResult q = par.align(p.ga, p.gb);
    EXPECT_EQ(q.snarl_alignments.size(), r.snarl_alignments.size());
    EXPECT_NEAR(q.total_score, r.total_score, 1e-3f);
    EXPECT_EQ(seqs(q.a_only_alleles), seqs(r.a_only_alleles));
    EXPECT_EQ(seqs(q.b_only_alleles), seqs(r.b_only_alleles));
}

TEST(Pairing, IntervalPairingWithoutLinkersIsOrderPreserving) {
    // No linkers.  The degree rule still admits the allele nodes of simple
    // bubbles (in- and out-degree 1), so shared alleles act as anchors and
    // every shared slot anchors itself; private slots fall into the preceding
    // interval and stay unpaired.  Either way the DP must recover the truth.
    Pair p = build("SASSBSSS", /*distinct=*/false, /*linker_every=*/0);
    GraphAligner seq(params(PairingMode::Interval));
    AlignmentResult r = seq.align(p.ga, p.gb);
    EXPECT_EQ(r.snarl_alignments.size(), static_cast<std::size_t>(p.shared_slots));
    EXPECT_EQ(r.unpaired_snarls_a, 1u);
    EXPECT_EQ(r.unpaired_snarls_b, 1u);
    EXPECT_EQ(seqs(r.a_only_alleles), ms(p.truth_a));
    EXPECT_EQ(seqs(r.b_only_alleles), ms(p.truth_b));
}

TEST(Pairing, ModesAgreeWhenKeysAreDistinct) {
    Pair p = build("SASSBSSS", /*distinct=*/true, /*linker_every=*/0);
    AlignmentResult rb = GraphAligner(params(PairingMode::Boundary)).align(p.ga, p.gb);
    AlignmentResult ri = GraphAligner(params(PairingMode::Interval)).align(p.ga, p.gb);
    EXPECT_EQ(rb.snarl_alignments.size(), ri.snarl_alignments.size());
    EXPECT_NEAR(rb.total_score, ri.total_score, 1e-3f);
    EXPECT_EQ(seqs(rb.a_only_alleles), ms(p.truth_a));
    EXPECT_EQ(seqs(ri.a_only_alleles), ms(p.truth_a));
}
