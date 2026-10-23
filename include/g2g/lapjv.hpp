#pragma once
/**
 * include/g2g/lapjv.hpp
 * =====================
 * Self-contained Jonker-Volgenant shortest-augmenting-path algorithm for
 * the dense linear assignment problem (LAPJV).
 *
 * Reference:
 *   R. Jonker & A. Volgenant, "A Shortest Augmenting Path Algorithm for
 *   Dense and Sparse Linear Assignment Problems," Computing 38, 325–340 (1987).
 *
 * API
 * ---
 *   int lapjv(int n, const float* cost, int* row_sol, int* col_sol,
 *             float* u = nullptr, float* v = nullptr);
 *
 *   cost     : n×n cost matrix, row-major, cost[i*n + j] = cost of pairing row i with col j.
 *   row_sol  : output - row_sol[i] = column assigned to row i.
 *   col_sol  : output - col_sol[j] = row    assigned to col j.
 *   u, v     : optional dual variables (size n each); may be nullptr.
 *   returns  : 0 on success.
 *
 * For G2G usage the caller builds the cost matrix as
 *   C[i][j] = 1.0f - seq_similarity(allele_a[i], allele_b[j])
 * and pads rectangular problems to square with gap_cost rows/columns.
 *
 * This header-only implementation is intentionally simple and readable
 * rather than maximally optimised.  SIMD / AVX paths can be layered in
 * Phase 4 if profiling shows the matching inner loop as a bottleneck
 * (Phase 2 data: bipartite matching ≥97 % of runtime at k ≥ 7).
 */

#include <algorithm>
#include <cassert>
#include <limits>
#include <vector>

namespace g2g {

/**
 * lapjv - square assignment, minimise sum of cost[row_sol[i], i].
 *
 * Implementation follows the three-phase structure of Jonker & Volgenant:
 *   Phase 1: column reduction   - initial cheap assignment.
 *   Phase 2: reduction transfer - repair incomplete assignment.
 *   Phase 3: augmenting row reduction (shortest-path augmentation).
 */
inline int lapjv(int n,
                 const float* __restrict__ cost,
                 int*   __restrict__ row_sol,
                 int*   __restrict__ col_sol,
                 float* __restrict__ u = nullptr,
                 float* __restrict__ v = nullptr)
{
    assert(n > 0);

    constexpr float kInf = std::numeric_limits<float>::max() / 2.f;

    // Dual variables (u = row potentials, v = col potentials)
    std::vector<float> u_buf, v_buf;
    if (!u) { u_buf.assign(n, 0.f); u = u_buf.data(); }
    if (!v) { v_buf.assign(n, 0.f); v = v_buf.data(); }

    std::vector<int> p(n + 1, 0);   // p[j]: row assigned to col j (1-indexed rows; 0 = unassigned)
    std::vector<int> way(n + 1, 0); // augmenting path predecessor

    // Jonker-Volgenant: iterate over rows, find shortest augmenting path via
    // Dijkstra-like scan over columns, then augment.
    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;

        std::vector<float> minval(n + 1, kInf);
        std::vector<bool>  used(n + 1, false);

        do {
            used[j0] = true;
            int    i0   = p[j0];
            float  delta = kInf;
            int    j1    = -1;

            for (int j = 1; j <= n; ++j) {
                if (used[j]) continue;
                // reduced cost: cost[i0-1, j-1] - u[i0-1] - v[j-1]
                float val = cost[(i0 - 1) * n + (j - 1)] - u[i0 - 1] - v[j - 1];
                if (val < minval[j]) {
                    minval[j] = val;
                    way[j]    = j0;
                }
                if (minval[j] < delta) {
                    delta = minval[j];
                    j1    = j;
                }
            }

            // Update potentials.  Column j (1..n) stores its potential at
            // v[j-1], matching the reduced-cost convention above (line uses
            // v[j-1]).  j == 0 is the sentinel start column and has no v entry,
            // so only its row potential u[p[0]-1] is updated.  Updating v[j]
            // here instead would both read/write one past the n-element buffer
            // (heap overflow at j == n) and shift every column potential.
            for (int j = 0; j <= n; ++j) {
                if (used[j]) {
                    u[p[j] - 1] += delta;
                    if (j >= 1) v[j - 1] -= delta;
                } else {
                    minval[j] -= delta;
                }
            }

            j0 = j1;
        } while (p[j0] != 0);

        // Augment: trace back the path and flip assignments
        do {
            int j1 = way[j0];
            p[j0]  = p[j1];
            j0     = j1;
        } while (j0 != 0);
    }

    // Extract solution (convert from 1-indexed internals to 0-indexed output)
    for (int j = 1; j <= n; ++j) {
        if (p[j] != 0) {
            row_sol[p[j] - 1] = j - 1;
            col_sol[j - 1]    = p[j] - 1;
        }
    }
    return 0;
}

} // namespace g2g
