#!/usr/bin/env bash
# ============================================================
# Experiment 3 — Throughput Scaling
#
# Runs G2G on graphs of increasing size:
#   1K, 5K, 10K, 50K, 100K, 500K, 1M snarls
#
# Compares against three baselines:
#   Baseline 1: minimap2 pairwise asm5 (linearize both graphs)
#   Baseline 2: PanGraph (Julia, bacterial graph-to-graph)
#   Baseline 3: sequential G2G (single thread)
#
# Records: wall time, peak memory, snarl pairs, score
#
# Usage:
#   conda activate g2g_eval
#   bash experiments/exp3_throughput/run.sh
#
# Environment variables:
#   G2G_DATA_ROOT  data directory [default: $HOME/g2g_data]
#   THREADS        threads for G2G [default: all cores]
#   SKIP_PANGRAPH  skip PanGraph baseline [default: 0]
# ============================================================
set -euo pipefail

G2G_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
RESULTS="$G2G_ROOT/results/phase5/exp3"
BUILD="$G2G_ROOT/build/release"
G2G_BIN="$BUILD/g2g"
G2G_SEQ_BIN="$BUILD/g2g"   # sequential mode via --sequential flag
GEN_SCRIPT="$G2G_ROOT/experiments/exp3_throughput/gen_synthetic_gfa.py"
THREADS="${THREADS:-$(nproc)}"
SKIP_PANGRAPH="${SKIP_PANGRAPH:-0}"
TIMESTAMP=$(date +%Y%m%dT%H%M%S)

GFA_DIR="$DATA_ROOT/synthetic_gfa"
mkdir -p "$RESULTS" "$GFA_DIR"

echo "[exp3] G2G root:   $G2G_ROOT"
echo "[exp3] Results:    $RESULTS"
echo "[exp3] Threads:    $THREADS"
echo "[exp3] Timestamp:  $TIMESTAMP"
echo "[exp3] minimap2:   $(minimap2 --version 2>&1 | head -1)"

if [[ ! -f "$G2G_BIN" ]]; then
    echo "[exp3] Building G2G..."
    cmake -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -G Ninja "$G2G_ROOT"
    cmake --build "$BUILD" -j"$(nproc)"
fi

# ─────────────────────────────────────────────────────────────
# Graph sizes to test
# ─────────────────────────────────────────────────────────────
# Format: "n_snarls:seq_len:k_max"
# seq_len_mean=12 for SNP-scale, 100 for SV-scale mix
SIZES=(
    "1000:12:12"
    "5000:12:20"
    "10000:12:20"
    "50000:12:30"
    "100000:12:50"
    "500000:12:50"
    "1000000:12:50"
)
SEEDS=(42 123 456)

CSV="$RESULTS/exp3_throughput_${TIMESTAMP}.csv"
echo "n_snarls,seq_len,k_max,seed,tool,threads,wall_ms,peak_kb,snarl_pairs,score,exit_code" \
    > "$CSV"

# ─────────────────────────────────────────────────────────────
# Helper: run one tool with timing
# ─────────────────────────────────────────────────────────────
run_timed() {
    local label="$1"; local cmd="$2"
    local time_file="$3"; local out_file="$4"

    local start_ns; start_ns=$(date +%s%N)
    /usr/bin/time -v -o "$time_file" bash -c "$cmd" > "$out_file" 2>&1
    local rc=$?
    local end_ns; end_ns=$(date +%s%N)
    local wall_ms=$(( (end_ns - start_ns) / 1000000 ))

    local peak_kb=0
    peak_kb=$(grep "Maximum resident" "$time_file" 2>/dev/null | \
              grep -o '[0-9]*' | tail -1 || echo 0)

    echo "${wall_ms}:${peak_kb}:${rc}"
}

# ─────────────────────────────────────────────────────────────
# Main loop
# ─────────────────────────────────────────────────────────────
for size_spec in "${SIZES[@]}"; do
    IFS=':' read -r n sl kmax <<< "$size_spec"

for seed in "${SEEDS[@]}"; do
    echo ""
    echo "══ n_snarls=$n  seq_len=$sl  k_max=$kmax  seed=$seed ══"

    # ── Generate GFA files ────────────────────────────────────
    GFA_A="$GFA_DIR/ga_n${n}_sl${sl}_km${kmax}_s${seed}.gfa"
    GFA_B="$GFA_DIR/gb_n${n}_sl${sl}_km${kmax}_s${seed}.gfa"

    if [[ ! -f "$GFA_A" ]] || [[ ! -s "$GFA_A" ]]; then
        echo "[exp3] Generating GFA (n=$n, sl=$sl, kmax=$kmax, seed=$seed)..."
        GEN_START=$(date +%s%N)
        python3 "$GEN_SCRIPT" \
            --n-snarls "$n" \
            --seq-len "$sl" \
            --seq-len-sv 100 \
            --k-max "$kmax" \
            --seed "$seed" \
            --out-a "$GFA_A" \
            --out-b "$GFA_B" \
            2>&1 | tee -a "$RESULTS/gen_log_${TIMESTAMP}.txt"
        GEN_MS=$(( ($(date +%s%N) - GEN_START) / 1000000 ))
        echo "[exp3] Generated in ${GEN_MS} ms"
    else
        echo "[exp3] GFA files exist, skipping generation"
    fi

    TMPDIR_RUN=$(mktemp -d "$RESULTS/run_${n}_${seed}_XXXXXX")
    trap "rm -rf $TMPDIR_RUN" EXIT

    # ── Baseline 1: minimap2 (linearize + align) ─────────────
    echo "[exp3] Baseline 1: minimap2..."

    # Extract haplotype FASTA from GFA
    awk '/^P/{
        name=$2; nodes=$3; split(nodes,a,",");
        seq=""; for(i in a){n=a[i]; sub(/[+-]/,"",n); seq=seq n}
        print ">"name"\n"seq
    }' "$GFA_A" | \
    awk 'NR==FNR{if(/^S/){s[$2]=$3} next}
         /^>/{print; next}
         {gsub(/[^ACGTN]/, "", $0); while(length($0)>0){print substr($0,1,80); $0=substr($0,81)} }' \
    /dev/null - 2>/dev/null \
    > "$TMPDIR_RUN/hap_a.fa" || true

    python3 - "$GFA_A" "$TMPDIR_RUN/hap_a.fa" <<'PYEOF'
import sys
seqs = {}
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('S'):
            p = line.split('\t')
            seqs[p[1]] = p[2].strip()
paths = {}
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('P'):
            p = line.split('\t')
            name = p[1]
            nodes = [n.rstrip('+-') for n in p[2].split(',')]
            paths[name] = ''.join(seqs.get(n,'') for n in nodes)
with open(sys.argv[2], 'w') as out:
    for name, seq in paths.items():
        out.write(f'>{name}\n{seq}\n')
PYEOF

    python3 - "$GFA_B" "$TMPDIR_RUN/hap_b.fa" <<'PYEOF'
import sys
seqs = {}
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('S'):
            p = line.split('\t')
            seqs[p[1]] = p[2].strip()
paths = {}
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('P'):
            p = line.split('\t')
            name = p[1]
            nodes = [n.rstrip('+-') for n in p[2].split(',')]
            paths[name] = ''.join(seqs.get(n,'') for n in nodes)
with open(sys.argv[2], 'w') as out:
    for name, seq in paths.items():
        out.write(f'>{name}\n{seq}\n')
PYEOF

    MM2_CMD="minimap2 -x asm5 --secondary=no -t $THREADS $TMPDIR_RUN/hap_b.fa $TMPDIR_RUN/hap_a.fa"
    result=$(run_timed "minimap2" "$MM2_CMD" \
        "$TMPDIR_RUN/time_mm2.txt" "$TMPDIR_RUN/mm2_out.paf")
    IFS=':' read -r mm2_wall mm2_mem mm2_rc <<< "$result"
    echo "  minimap2: ${mm2_wall}ms  ${mm2_mem}KB  rc=$mm2_rc"
    echo "$n,$sl,$kmax,$seed,minimap2,$THREADS,$mm2_wall,$mm2_mem,0,0,$mm2_rc" >> "$CSV"

    # ── Baseline 2: PanGraph (skip if Julia not available) ───
    if [[ "$SKIP_PANGRAPH" == "0" ]] && command -v pangraph &>/dev/null; then
        echo "[exp3] Baseline 2: PanGraph..."
        # Convert GFA to FASTA (PanGraph works on FASTA genomes)
        awk '/^P/{
            name=$2; nodes=$3; split(nodes,a,",");
            printf ">"name"\n"
            for(i=1;i<=length(a);i++){n=a[i]; sub(/[+-]/,"",n); printf smap[n]}
            printf "\n"
        }
        /^S/{smap[$2]=$3}' "$GFA_A" > "$TMPDIR_RUN/genome_a.fa" 2>/dev/null || true

        awk '/^P/{
            name=$2; nodes=$3; split(nodes,a,",");
            printf ">"name"\n"
            for(i=1;i<=length(a);i++){n=a[i]; sub(/[+-]/,"",n); printf smap[n]}
            printf "\n"
        }
        /^S/{smap[$2]=$3}' "$GFA_B" > "$TMPDIR_RUN/genome_b.fa" 2>/dev/null || true

        PG_CMD="pangraph build --circular false \
            $TMPDIR_RUN/genome_a.fa $TMPDIR_RUN/genome_b.fa \
            -o $TMPDIR_RUN/pangraph_out.json"
        result=$(run_timed "pangraph" "$PG_CMD" \
            "$TMPDIR_RUN/time_pg.txt" "$TMPDIR_RUN/pg_log.txt" || echo "0:0:1")
        IFS=':' read -r pg_wall pg_mem pg_rc <<< "$result"
        echo "  PanGraph: ${pg_wall}ms  ${pg_mem}KB  rc=$pg_rc"
        echo "$n,$sl,$kmax,$seed,pangraph,1,$pg_wall,$pg_mem,0,0,$pg_rc" >> "$CSV"
    else
        echo "[exp3] Skipping PanGraph (not installed or SKIP_PANGRAPH=1)"
        echo "$n,$sl,$kmax,$seed,pangraph,0,0,0,0,0,-1" >> "$CSV"
    fi

    # ── Baseline 3: Sequential G2G ────────────────────────────
    echo "[exp3] Baseline 3: Sequential G2G..."
    SEQ_CMD="$G2G_BIN $GFA_A $GFA_B --sequential \
        --gaf $TMPDIR_RUN/seq_out.gaf --vcf $TMPDIR_RUN/seq_out.vcf"
    result=$(run_timed "g2g_seq" "$SEQ_CMD" \
        "$TMPDIR_RUN/time_seq.txt" "$TMPDIR_RUN/g2g_seq_log.txt")
    IFS=':' read -r seq_wall seq_mem seq_rc <<< "$result"

    seq_pairs=$(grep -c "^" "$TMPDIR_RUN/seq_out.gaf" 2>/dev/null || echo 0)
    seq_score=$(grep "total_score=" "$TMPDIR_RUN/g2g_seq_log.txt" 2>/dev/null | \
                grep -o 'total_score=[0-9.]*' | head -1 | cut -d= -f2 || echo 0)
    echo "  G2G-seq: ${seq_wall}ms  ${seq_mem}KB  pairs=$seq_pairs  score=$seq_score"
    echo "$n,$sl,$kmax,$seed,g2g_sequential,1,$seq_wall,$seq_mem,$seq_pairs,$seq_score,$seq_rc" >> "$CSV"

    # ── Production G2G (parallel) ─────────────────────────────
    echo "[exp3] Production G2G (parallel, $THREADS threads)..."
    PAR_CMD="$G2G_BIN $GFA_A $GFA_B --threads $THREADS \
        --gaf $TMPDIR_RUN/par_out.gaf --vcf $TMPDIR_RUN/par_out.vcf"
    result=$(run_timed "g2g_par" "$PAR_CMD" \
        "$TMPDIR_RUN/time_par.txt" "$TMPDIR_RUN/g2g_par_log.txt")
    IFS=':' read -r par_wall par_mem par_rc <<< "$result"

    par_pairs=$(grep -c "^" "$TMPDIR_RUN/par_out.gaf" 2>/dev/null || echo 0)
    par_score=$(grep "total_score=" "$TMPDIR_RUN/g2g_par_log.txt" 2>/dev/null | \
                grep -o 'total_score=[0-9.]*' | head -1 | cut -d= -f2 || echo 0)
    echo "  G2G-par: ${par_wall}ms  ${par_mem}KB  pairs=$par_pairs  score=$par_score"
    echo "$n,$sl,$kmax,$seed,g2g_parallel,$THREADS,$par_wall,$par_mem,$par_pairs,$par_score,$par_rc" >> "$CSV"

    rm -rf "$TMPDIR_RUN"

done; done

echo ""
echo "[exp3] Results written to: $CSV"
echo "[exp3] Summary:"
python3 - "$CSV" <<'PYEOF'
import csv, sys
with open(sys.argv[1]) as f:
    rows = list(csv.DictReader(f))
tools = sorted(set(r['tool'] for r in rows))
sizes = sorted(set(int(r['n_snarls']) for r in rows))
print(f"{'n_snarls':>12}", end="")
for t in tools: print(f"  {t:>18}", end="")
print()
for n in sizes:
    print(f"{n:>12}", end="")
    for t in tools:
        subset = [r for r in rows if int(r['n_snarls'])==n and r['tool']==t]
        if subset:
            avg = sum(int(r['wall_ms']) for r in subset) / len(subset)
            print(f"  {avg:>16.0f}ms", end="")
        else:
            print(f"  {'N/A':>18}", end="")
    print()
PYEOF
