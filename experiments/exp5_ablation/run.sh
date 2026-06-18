#!/usr/bin/env bash
# ============================================================
# Experiment 5 — Ablation Study
#
# Measures the isolated contribution of each design decision:
#
#   A1: Cache-oblivious vEB layout vs. BFS layout
#       → L3 miss rate, throughput
#   A2: Work-stealing vs. static partition
#       → load imbalance, throughput
#   A3: Steal-largest vs. steal-random
#       → throughput on heavy-tailed workload
#   A4: Minimizer index vs. naïve O(|VA||VB|) anchor detection
#       → anchor detection time, throughput
#   A5: Hungarian matching vs. brute-force (at k=5,10,20,50)
#       → matching time per snarl pair
#
# For A1–A3: graph size = 10K snarls (representative of HPRC)
# For A4:    graph size = 100K snarls (reveals index benefits)
# For A5:    synthetic micro-benchmarks with fixed k
#
# Usage:
#   conda activate g2g_eval
#   bash experiments/exp5_ablation/run.sh
#
# Output: results/phase5/exp5/
# ============================================================
set -euo pipefail

G2G_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
RESULTS="$G2G_ROOT/results/phase5/exp5"
BUILDS="$G2G_ROOT/build"
GEN_SCRIPT="$G2G_ROOT/experiments/exp3_throughput/gen_synthetic_gfa.py"
THREADS="${THREADS:-$(nproc)}"
REPEATS="${REPEATS:-5}"
TIMESTAMP=$(date +%Y%m%dT%H%M%S)

GFA_DIR="$DATA_ROOT/synthetic_gfa"
mkdir -p "$RESULTS" "$GFA_DIR"

echo "[exp5] Results:  $RESULTS"
echo "[exp5] Threads:  $THREADS"
echo "[exp5] Repeats:  $REPEATS"

# ─────────────────────────────────────────────────────────────
# Verify all ablation binaries are built
# ─────────────────────────────────────────────────────────────
VARIANTS=(
    "release"
    "abl_bfs_layout"
    "abl_steal_random"
    "abl_naive_anchor"
    "abl_brute_matching"
    "abl_static_partition"
)

for v in "${VARIANTS[@]}"; do
    bin="$BUILDS/$v/g2g"
    if [[ ! -f "$bin" ]]; then
        echo "[exp5] Missing binary: $bin — run 01_build_all.sh first"
        exit 1
    fi
done
echo "[exp5] All ${#VARIANTS[@]} ablation binaries found."

# ─────────────────────────────────────────────────────────────
# Helper: run g2g with timing, return "wall_ms:peak_kb:rc"
# ─────────────────────────────────────────────────────────────
run_g2g() {
    local bin="$1"; local gfa_a="$2"; local gfa_b="$3"
    local extra_flags="$4"
    local time_file="$5"; local out_gaf="$6"
    local start end

    start=$(date +%s%N)
    /usr/bin/time -v -o "$time_file" \
        "$bin" "$gfa_a" "$gfa_b" $extra_flags \
        --gaf "$out_gaf" \
        > "${time_file%.txt}_log.txt" 2>&1
    local rc=$?
    end=$(date +%s%N)
    local wall_ms=$(( (end - start) / 1000000 ))
    local peak_kb
    peak_kb=$(grep "Maximum resident" "$time_file" 2>/dev/null | \
              grep -o '[0-9]*' | tail -1 || echo 0)

    echo "${wall_ms}:${peak_kb}:${rc}"
}

# ─────────────────────────────────────────────────────────────
# A1–A3: 10K snarl graphs
# ─────────────────────────────────────────────────────────────
N_MED=10000
SEEDS_MED=(42 123 456 789 1234)

CSV_MAIN="$RESULTS/exp5_ablation_main_${TIMESTAMP}.csv"
echo "ablation,variant,n_snarls,seed,threads,repeat,wall_ms,peak_kb,rc" \
    > "$CSV_MAIN"

for seed in "${SEEDS_MED[@]}"; do
    GFA_A="$GFA_DIR/ga_n${N_MED}_s${seed}.gfa"
    GFA_B="$GFA_DIR/gb_n${N_MED}_s${seed}.gfa"
    if [[ ! -f "$GFA_A" ]]; then
        echo "[exp5] Generating n=$N_MED seed=$seed..."
        python3 "$GEN_SCRIPT" \
            --n-snarls "$N_MED" --seq-len 12 --seq-len-sv 100 \
            --k-max 50 --seed "$seed" \
            --out-a "$GFA_A" --out-b "$GFA_B"
    fi
done

# A1: Cache-oblivious vEB layout vs. BFS layout
echo "[exp5] A1: vEB vs BFS layout..."
for seed in "${SEEDS_MED[@]}"; do
    GFA_A="$GFA_DIR/ga_n${N_MED}_s${seed}.gfa"
    GFA_B="$GFA_DIR/gb_n${N_MED}_s${seed}.gfa"
    for rep in $(seq 1 "$REPEATS"); do
        TMPD=$(mktemp -d "$RESULTS/abl_XXXXXX")
        trap "rm -rf $TMPD" EXIT

        result=$(run_g2g "$BUILDS/release/g2g" \
            "$GFA_A" "$GFA_B" "--threads $THREADS" \
            "$TMPD/time_veb.txt" "$TMPD/out_veb.gaf")
        IFS=':' read -r w m rc <<< "$result"
        echo "A1,veb_layout,$N_MED,$seed,$THREADS,$rep,$w,$m,$rc" >> "$CSV_MAIN"

        result=$(run_g2g "$BUILDS/abl_bfs_layout/g2g" \
            "$GFA_A" "$GFA_B" "--threads $THREADS" \
            "$TMPD/time_bfs.txt" "$TMPD/out_bfs.gaf")
        IFS=':' read -r w m rc <<< "$result"
        echo "A1,bfs_layout,$N_MED,$seed,$THREADS,$rep,$w,$m,$rc" >> "$CSV_MAIN"

        rm -rf "$TMPD"
    done
done
echo "[exp5] A1 done."

# A2: Work-stealing vs. static partition
echo "[exp5] A2: Work-stealing vs static partition..."
for seed in "${SEEDS_MED[@]}"; do
    GFA_A="$GFA_DIR/ga_n${N_MED}_s${seed}.gfa"
    GFA_B="$GFA_DIR/gb_n${N_MED}_s${seed}.gfa"
    for rep in $(seq 1 "$REPEATS"); do
        TMPD=$(mktemp -d "$RESULTS/abl_XXXXXX")
        trap "rm -rf $TMPD" EXIT

        result=$(run_g2g "$BUILDS/release/g2g" \
            "$GFA_A" "$GFA_B" "--threads $THREADS" \
            "$TMPD/time_ws.txt" "$TMPD/out_ws.gaf")
        IFS=':' read -r w m rc <<< "$result"
        echo "A2,work_stealing,$N_MED,$seed,$THREADS,$rep,$w,$m,$rc" >> "$CSV_MAIN"

        result=$(run_g2g "$BUILDS/abl_static_partition/g2g" \
            "$GFA_A" "$GFA_B" "--threads $THREADS" \
            "$TMPD/time_sp.txt" "$TMPD/out_sp.gaf")
        IFS=':' read -r w m rc <<< "$result"
        echo "A2,static_partition,$N_MED,$seed,$THREADS,$rep,$w,$m,$rc" >> "$CSV_MAIN"

        rm -rf "$TMPD"
    done
done
echo "[exp5] A2 done."

# A3: Steal-largest vs. steal-random
echo "[exp5] A3: Steal-largest vs steal-random..."
for seed in "${SEEDS_MED[@]}"; do
    GFA_A="$GFA_DIR/ga_n${N_MED}_s${seed}.gfa"
    GFA_B="$GFA_DIR/gb_n${N_MED}_s${seed}.gfa"
    for rep in $(seq 1 "$REPEATS"); do
        TMPD=$(mktemp -d "$RESULTS/abl_XXXXXX")
        trap "rm -rf $TMPD" EXIT

        result=$(run_g2g "$BUILDS/release/g2g" \
            "$GFA_A" "$GFA_B" "--threads $THREADS" \
            "$TMPD/time_sl.txt" "$TMPD/out_sl.gaf")
        IFS=':' read -r w m rc <<< "$result"
        echo "A3,steal_largest,$N_MED,$seed,$THREADS,$rep,$w,$m,$rc" >> "$CSV_MAIN"

        result=$(run_g2g "$BUILDS/abl_steal_random/g2g" \
            "$GFA_A" "$GFA_B" "--threads $THREADS" \
            "$TMPD/time_sr.txt" "$TMPD/out_sr.gaf")
        IFS=':' read -r w m rc <<< "$result"
        echo "A3,steal_random,$N_MED,$seed,$THREADS,$rep,$w,$m,$rc" >> "$CSV_MAIN"

        rm -rf "$TMPD"
    done
done
echo "[exp5] A3 done."

# ─────────────────────────────────────────────────────────────
# A4: Minimizer index vs. naïve anchor detection (100K snarls)
# ─────────────────────────────────────────────────────────────
N_LARGE=100000
CSV_A4="$RESULTS/exp5_ablation_anchor_${TIMESTAMP}.csv"
echo "ablation,variant,n_snarls,seed,threads,repeat,wall_ms,peak_kb,rc" \
    > "$CSV_A4"

for seed in "${SEEDS_MED[@]:0:3}"; do  # 3 seeds for large graphs
    GFA_A="$GFA_DIR/ga_n${N_LARGE}_s${seed}.gfa"
    GFA_B="$GFA_DIR/gb_n${N_LARGE}_s${seed}.gfa"
    if [[ ! -f "$GFA_A" ]]; then
        echo "[exp5] Generating n=$N_LARGE seed=$seed..."
        python3 "$GEN_SCRIPT" \
            --n-snarls "$N_LARGE" --seq-len 12 --seq-len-sv 100 \
            --k-max 50 --seed "$seed" \
            --out-a "$GFA_A" --out-b "$GFA_B"
    fi

    echo "[exp5] A4: n=$N_LARGE seed=$seed..."
    for rep in $(seq 1 3); do  # fewer repeats for large graphs
        TMPD=$(mktemp -d "$RESULTS/abl_XXXXXX")
        trap "rm -rf $TMPD" EXIT

        result=$(run_g2g "$BUILDS/release/g2g" \
            "$GFA_A" "$GFA_B" "--threads $THREADS" \
            "$TMPD/time_mi.txt" "$TMPD/out_mi.gaf")
        IFS=':' read -r w m rc <<< "$result"
        echo "A4,minimizer_index,$N_LARGE,$seed,$THREADS,$rep,$w,$m,$rc" >> "$CSV_A4"

        result=$(run_g2g "$BUILDS/abl_naive_anchor/g2g" \
            "$GFA_A" "$GFA_B" "--threads $THREADS" \
            "$TMPD/time_na.txt" "$TMPD/out_na.gaf")
        IFS=':' read -r w m rc <<< "$result"
        echo "A4,naive_anchor,$N_LARGE,$seed,$THREADS,$rep,$w,$m,$rc" >> "$CSV_A4"

        rm -rf "$TMPD"
    done
done
echo "[exp5] A4 done."

# ─────────────────────────────────────────────────────────────
# A5: Hungarian vs. brute-force matching at fixed k
# ─────────────────────────────────────────────────────────────
CSV_A5="$RESULTS/exp5_ablation_matching_${TIMESTAMP}.csv"
echo "ablation,variant,k,n_snarls,seed,threads,repeat,wall_ms,peak_kb,rc" \
    > "$CSV_A5"

K_VALUES=(5 6 7 8 10 15 20)
N_A5=1000  # small graph, only a few snarls — isolate matching cost

for k in "${K_VALUES[@]}"; do
    for seed in 42 123 456; do
        # Generate a graph where ALL snarls have exactly k alleles each
        GFA_A="$GFA_DIR/ga_k${k}_n${N_A5}_s${seed}.gfa"
        GFA_B="$GFA_DIR/gb_k${k}_n${N_A5}_s${seed}.gfa"
        if [[ ! -f "$GFA_A" ]]; then
            python3 - "$k" "$N_A5" "$seed" "$GFA_A" "$GFA_B" <<'PYEOF'
import sys, random
k, n_snarls, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
out_a, out_b = sys.argv[4], sys.argv[5]
rng = random.Random(seed)
DNA = "ACGT"
def rseq(l): return ''.join(rng.choice(DNA) for _ in range(l))

nodes_a, nodes_b = {}, {}
edges_a, edges_b = [], []
path_a, path_b = [], []

for i in range(n_snarls + 1):
    s = rseq(8)
    n = f"anc_{i}"
    nodes_a[n] = nodes_b[n] = s
path_a.append("anc_0"); path_b.append("anc_0")

for s in range(n_snarls):
    src, snk = f"anc_{s}", f"anc_{s+1}"
    for j in range(k):
        seq = rseq(12)
        na, nb = f"s{s}_a{j}", f"s{s}_b{j}"
        nodes_a[na] = seq; nodes_b[nb] = seq
        edges_a += [(src,na),(na,snk)]
        edges_b += [(src,nb),(nb,snk)]
    path_a += [f"s{s}_a0", snk]
    path_b += [f"s{s}_b0", snk]

for out_f, nodes, edges, path in [
    (out_a, nodes_a, edges_a, path_a),
    (out_b, nodes_b, edges_b, path_b),
]:
    with open(out_f, 'w') as f:
        f.write("H\tVN:Z:1.1\n")
        for n,seq in nodes.items(): f.write(f"S\t{n}\t{seq}\n")
        seen = set()
        for u,v in edges:
            if (u,v) not in seen:
                f.write(f"L\t{u}\t+\t{v}\t+\t0M\n"); seen.add((u,v))
        segs = ",".join(f"{n}+" for n in path)
        ovlps = ",".join("0M" for _ in range(len(path)-1))
        f.write(f"P\thap\t{segs}\t{ovlps}\n")
PYEOF
        fi

        echo "[exp5] A5: k=$k  seed=$seed..."
        for rep in $(seq 1 "$REPEATS"); do
            TMPD=$(mktemp -d "$RESULTS/abl_XXXXXX")
            trap "rm -rf $TMPD" EXIT

            result=$(run_g2g "$BUILDS/release/g2g" \
                "$GFA_A" "$GFA_B" "--sequential" \
                "$TMPD/time_hung.txt" "$TMPD/out_hung.gaf")
            IFS=':' read -r w m rc <<< "$result"
            echo "A5,hungarian,$k,$N_A5,$seed,1,$rep,$w,$m,$rc" >> "$CSV_A5"

            if [[ "$k" -le 8 ]]; then  # brute-force only feasible for k<=8
                result=$(run_g2g "$BUILDS/abl_brute_matching/g2g" \
                    "$GFA_A" "$GFA_B" "--sequential" \
                    "$TMPD/time_brute.txt" "$TMPD/out_brute.gaf")
                IFS=':' read -r w m rc <<< "$result"
                echo "A5,brute_force,$k,$N_A5,$seed,1,$rep,$w,$m,$rc" >> "$CSV_A5"
            else
                echo "A5,brute_force,$k,$N_A5,$seed,1,$rep,timeout,0,skip" >> "$CSV_A5"
            fi
            rm -rf "$TMPD"
        done
    done
done
echo "[exp5] A5 done."

# ─────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────
echo ""
echo "[exp5] Results:"
echo "  $CSV_MAIN"
echo "  $CSV_A4"
echo "  $CSV_A5"
echo ""
echo "[exp5] Ablation speedup summary (A1–A3, n=$N_MED):"
python3 - "$CSV_MAIN" <<'PYEOF'
import csv, sys, collections
rows = list(csv.DictReader(open(sys.argv[1])))
for ab in ['A1','A2','A3']:
    subset = [r for r in rows if r['ablation']==ab]
    variants = sorted(set(r['variant'] for r in subset))
    print(f"\n  {ab}:")
    for v in variants:
        ws = [int(r['wall_ms']) for r in subset if r['variant']==v and r['wall_ms'].isdigit()]
        if ws:
            avg = sum(ws)/len(ws)
            print(f"    {v:25s}: {avg:8.1f} ms (n={len(ws)})")
PYEOF
