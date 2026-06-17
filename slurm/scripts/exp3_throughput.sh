#!/bin/bash --login
# =============================================================================
# Experiment 3 — Throughput Scaling & Tool Comparison   (Q2 + Q3)
# =============================================================================
#SBATCH --job-name=g2g_exp3
#SBATCH --output=slurm/logs/exp3-%j.out
#SBATCH --error=slurm/logs/exp3-%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
#SBATCH --constraint=cpu_amd_epyc_7763

set -euo pipefail

# ── Paths (override via environment) ─────────────────────────────────────────
if [[ -n "${G2G_REPO:-}" ]]; then
    REPO="$G2G_REPO"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -e "$SLURM_SUBMIT_DIR/slurm/scripts/exp3_throughput.sh" ]]; then
    REPO="$SLURM_SUBMIT_DIR"
else
    _src="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO="$(cd "$_src/../.." && pwd)"
fi
if [[ ! -e "$REPO/slurm/scripts/exp3_throughput.sh" ]]; then
    echo "FATAL: could not locate the G2G repo (REPO='$REPO')." >&2
    echo "  FIX: sbatch --export=ALL,G2G_REPO=\$PWD slurm/scripts/exp3_throughput.sh" >&2
    exit 1
fi
SCRIPT_DIR="$REPO/slurm/scripts"
BUILD="$REPO/build"
G2G="$BUILD/g2g"
GEN="$REPO/slurm/synth/gen_hprc_scale_gfa.py"
LINEARISE="$REPO/slurm/synth/linearise_gfa.py"
OUTDIR="$REPO/results/exp3"
GRAPHDIR="${G2G_GRAPHDIR:-$REPO/results/exp3/graphs}"
RUNS="${G2G_RUNS:-3}"

# Graph-size sweep. PARALLEL runs at every size here.
SIZES_DEFAULT="1000 10000 100000 1000000 10000000"
SIZES="${G2G_SIZES:-$SIZES_DEFAULT}"

# ── SEQUENTIAL-BASELINE SIZE GATE (the runaway fix) ──────────────────────────
# The single-threaded baseline is O(n_A * n_B) — QUADRATIC in snarl count.
# Measured: n=10k -> ~4815 s/run. Extrapolated: n=10^5 -> ~5.6 days/run,
# n=10^6 -> ~1.5 years/run. Running it at large n is BOTH infeasible and
# unnecessary: it is only needed as the speedup *denominator* at sizes where it
# is tractable. Above this threshold we run PARALLEL ONLY and report throughput
# scaling (the headline figure), not a sequential-vs-parallel ratio.
#   Default 10000 keeps each sequential run ~80 min and the whole sequential
#   set bounded. Override with G2G_SEQ_MAX_N if you have time/headroom.
SEQ_MAX_N="${G2G_SEQ_MAX_N:-10000}"

# Haplotype-panel sweep (Exp 3b): fixed graph size, growing panel.
HAP_N="${G2G_HAP_N:-10000}"
HAP_PANELS="${G2G_HAP_PANELS:-2 10 30 100 300}"

BASELINE_TIMEOUT="${G2G_BASELINE_TIMEOUT:-7200}"   # 2 h

NCPU="$(nproc)"

mkdir -p "$OUTDIR" "$GRAPHDIR" "$REPO/slurm/logs"
source "$SCRIPT_DIR/lib_measure.sh"

if [[ -n "${CONDA_BASE:-}" && -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
fi
conda activate "${G2G_ENV:-phase5}" 2>/dev/null || log "WARN: could not activate conda env"

capture_env "$OUTDIR/env.txt"
log "exp3 start  repo=$REPO  binary=$G2G  sizes=[$SIZES]  runs=$RUNS  seq_max_n=$SEQ_MAX_N"

# ── Rebuild Release ──────────────────────────────────────────────────────────
if command -v cmake >/dev/null 2>&1; then
    if cmake -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_CXX_COMPILER="$(which g++)" \
          -DCMAKE_C_COMPILER="$(which gcc)" \
          -DCMAKE_CXX_FLAGS="-O3 -march=native" >/dev/null 2>&1 \
       && cmake --build "$BUILD" -j"$NCPU" --target g2g >/dev/null 2>&1; then
        log "binary rebuilt (Release, -O3 -march=native)"
    else
        log "WARN: rebuild failed; using existing binary at $G2G"
    fi
else
    log "WARN: cmake not found; using existing binary at $G2G"
fi
if [[ ! -x "$G2G" ]]; then
    log "FATAL: no g2g binary at $G2G — run setup_cluster.sh first"
    exit 1
fi

VERBOSE_FLAG=""
if "$G2G" --help 2>&1 | grep -q -- "--verbose-timing"; then
    VERBOSE_FLAG="--verbose-timing"
    log "binary supports --verbose-timing: per-phase breakdown ENABLED"
else
    log "binary lacks --verbose-timing: per-phase columns will be NA"
fi

have() { command -v "$1" >/dev/null 2>&1; }
HAVE_VG=0;          have vg          && HAVE_VG=1
HAVE_GRAPHALIGNER=0; have GraphAligner && HAVE_GRAPHALIGNER=1
HAVE_MINIGRAPH=0;   have minigraph   && HAVE_MINIGRAPH=1
HAVE_MINIMAP2=0;    have minimap2    && HAVE_MINIMAP2=1
log "baselines: vg=$HAVE_VG GraphAligner=$HAVE_GRAPHALIGNER minigraph=$HAVE_MINIGRAPH minimap2=$HAVE_MINIMAP2"

parse_phases() {
    local logfile="$1"
    local line
    line=$(grep -m1 "\[timing\]" "$logfile" 2>/dev/null || true)
    if [[ -z "$line" ]]; then echo "NA NA NA NA"; return; fi
    local io idx dec dp
    io=$( echo "$line" | grep -oP 'io=\K[0-9.]+'        || echo NA)
    idx=$(echo "$line" | grep -oP 'index=\K[0-9.]+'     || echo NA)
    dec=$(echo "$line" | grep -oP 'decompose=\K[0-9.]+' || echo NA)
    dp=$( echo "$line" | grep -oP 'dp=\K[0-9.]+'        || echo NA)
    echo "${io:-NA} ${idx:-NA} ${dec:-NA} ${dp:-NA}"
}

RESULTS="$OUTDIR/throughput_scaling.tsv"
FAILMODE="$OUTDIR/baseline_failure_modes.tsv"
MANIFEST="$OUTDIR/manifest.sha256"
: > "$MANIFEST"
tsv_init "$RESULTS" "n_snarls\tsystem\trun\twall_sec\tpeak_rss_kb\tuser_sec\tsys_sec\tcpu_pct\texit\tphase_io\tphase_index\tphase_decompose\tphase_dp"
tsv_init "$FAILMODE" "n_snarls\ttool\tfailure_mode\tdetail"

run_g2g() {
    local n="$1" label="$2" run="$3"; shift 3
    local ga="$GRAPHDIR/ga_${n}.gfa" gb="$GRAPHDIR/gb_${n}.gfa"
    local g2glog; g2glog=$(mktemp)
    local metrics
    metrics=$(measure_full_capture "$g2glog" \
              "$G2G" "$ga" "$gb" --gaf /dev/null $VERBOSE_FLAG "$@")
    read -r wall peak user sys cpu ec <<< "$metrics"
    local phases; phases=$(parse_phases "$g2glog")
    read -r p_io p_idx p_dec p_dp <<< "$phases"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$n" "$label" "$run" "$wall" "$peak" "$user" "$sys" "$cpu" "$ec" \
        "$p_io" "$p_idx" "$p_dec" "$p_dp" >> "$RESULTS"
    log "  [$label] n=$n run=$run wall=${wall}s rss=${peak}kB cpu=${cpu}% exit=$ec"
    rm -f "$g2glog"
    return 0
}

record_failure() {
    local n="$1" tool="$2" mode="$3" detail="$4"
    printf '%s\t%s\t%s\t%s\n' "$n" "$tool" "$mode" "$detail" >> "$FAILMODE"
    log "  [$tool] n=$n INFEASIBLE: $mode ($detail)"
}

linearise() {
    local gfa="$1" fasta="$2"
    if [[ -f "$LINEARISE" ]]; then
        python3 "$LINEARISE" "$gfa" "$fasta"
    else
        python3 - "$gfa" "$fasta" <<'PYEOF'
import sys
def rc(s):
    t={'A':'T','T':'A','C':'G','G':'C','a':'t','t':'a','c':'g','g':'c'}
    return ''.join(t.get(c,'N') for c in reversed(s))
gfa,out=sys.argv[1],sys.argv[2]
seqs={}; paths={}
with open(gfa) as f:
    for line in f:
        if line[0]=='S':
            p=line.rstrip('\n').split('\t'); seqs[p[1]]=p[2] if p[2]!='*' else ''
        elif line[0]=='P':
            p=line.rstrip('\n').split('\t'); name=p[1]
            segs=p[2].split(','); s=''
            for seg in segs:
                o=seg[-1]; nm=seg[:-1] if o in '+-' else seg
                bp=seqs.get(nm,'')
                s+=bp if o!='-' else rc(bp)
            paths[name]=s
with open(out,'w') as f:
    for name,s in paths.items():
        f.write(f'>{name}\n{s}\n')
PYEOF
    fi
}

run_baseline() {
    local n="$1" tool="$2"; shift 2
    local blog; blog=$(mktemp)
    local metrics
    metrics=$(measure_full timeout "$BASELINE_TIMEOUT" "$@" 2>"$blog")
    read -r wall peak user sys cpu ec <<< "$metrics"
    if [[ "$ec" == "124" ]]; then
        record_failure "$n" "$tool" "timeout" ">${BASELINE_TIMEOUT}s wall"
    elif [[ "$ec" != "0" ]]; then
        if [[ "$ec" == "137" ]] || grep -qiE "out of memory|bad_alloc|cannot allocate" "$blog"; then
            record_failure "$n" "$tool" "oom" "exit=$ec peak_rss=${peak}kB"
        else
            local snippet; snippet=$(tail -1 "$blog" | tr '\t' ' ' | cut -c1-80)
            record_failure "$n" "$tool" "error" "exit=$ec '$snippet'"
        fi
    else
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$n" "$tool" "1" "$wall" "$peak" "$user" "$sys" "$cpu" "$ec" \
            "NA" "NA" "NA" "NA" >> "$RESULTS"
        log "  [$tool] n=$n wall=${wall}s rss=${peak}kB cpu=${cpu}%"
    fi
    rm -f "$blog"
    return 0
}

# =============================================================================
# PART 1 — THROUGHPUT SWEEP
# =============================================================================
log "=== PART 1: throughput sweep ==="
for N in $SIZES; do
    log "----- n_snarls=$N -----"
    GA="$GRAPHDIR/ga_${N}.gfa"; GB="$GRAPHDIR/gb_${N}.gfa"

    if [[ ! -s "$GA" || ! -s "$GB" ]]; then
        log "  generating graphs (n=$N)..."
        python3 "$GEN" --mode size --n_snarls "$N" --out_dir "$GRAPHDIR" --seed 42
    else
        log "  graphs present, reusing"
    fi
    echo "$(sha_of "$GA")  $(basename "$GA")" >> "$MANIFEST"
    echo "$(sha_of "$GB")  $(basename "$GB")" >> "$MANIFEST"

    for RUN in $(seq 1 "$RUNS"); do
        # ── Sequential baseline: ONLY at sizes where it is tractable ─────────
        if (( N <= SEQ_MAX_N )); then
            run_g2g "$N" "g2g_sequential" "$RUN" --sequential
        elif (( RUN == 1 )); then
            log "  [g2g_sequential] n=$N SKIPPED (n > G2G_SEQ_MAX_N=$SEQ_MAX_N; O(n^2) baseline infeasible — reporting parallel-only)"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$N" "g2g_sequential" "$RUN" "NA" "NA" "NA" "NA" "NA" "skipped_seq_max_n" \
                "NA" "NA" "NA" "NA" >> "$RESULTS"
        fi

        # ── Parallel: at EVERY size (this is the headline scaling curve) ─────
        run_g2g "$N" "g2g_parallel_${NCPU}t" "$RUN" --threads "$NCPU"
    done

    # ── Baselines ────────────────────────────────────────────────────────────
    FA="$GRAPHDIR/seqs_a_${N}.fa"; FB="$GRAPHDIR/seqs_b_${N}.fa"
    NEED_LINEARISE=0
    if [[ $HAVE_MINIMAP2 -eq 1 || $HAVE_GRAPHALIGNER -eq 1 || $HAVE_VG -eq 1 || $HAVE_MINIGRAPH -eq 1 ]]; then
        NEED_LINEARISE=1
    fi
    if [[ $NEED_LINEARISE -eq 1 && ( ! -s "$FA" || ! -s "$FB" ) ]]; then
        log "  linearising graphs for baselines..."
        if ! timeout "$BASELINE_TIMEOUT" bash -c "$(declare -f linearise); linearise '$GA' '$FA'" 2>/dev/null; then
            record_failure "$N" "linearise" "path_explosion" "haplotype extraction exceeded ${BASELINE_TIMEOUT}s or OOM"
            FA=""; FB=""
        else
            linearise "$GB" "$FB" || { FA=""; FB=""; }
        fi
    fi

    if [[ -n "$FA" && -s "$FA" ]]; then
        [[ $HAVE_MINIMAP2 -eq 1 ]] && \
            run_baseline "$N" "minimap2_asm5" \
                minimap2 -x asm5 -t "$NCPU" --secondary=no "$FA" "$FB" -o /dev/null
        [[ $HAVE_GRAPHALIGNER -eq 1 ]] && \
            run_baseline "$N" "graphaligner" \
                GraphAligner -g "$GB" -f "$FA" -a /dev/null -x vg -t "$NCPU"
        if [[ $HAVE_VG -eq 1 ]]; then
            run_baseline "$N" "vg_giraffe" bash -c "
                set -e
                vg convert -g '$GB' -p > '$GRAPHDIR/gb_${N}.vg' 2>/dev/null
                vg index -x '$GRAPHDIR/gb_${N}.xg' '$GRAPHDIR/gb_${N}.vg' 2>/dev/null
                vg giraffe -x '$GRAPHDIR/gb_${N}.xg' -f '$FA' -t $NCPU > /dev/null 2>&1
            "
        fi
        [[ $HAVE_MINIGRAPH -eq 1 ]] && \
            run_baseline "$N" "minigraph" \
                minigraph -cxggs -t "$NCPU" "$GB" "$FA"
    else
        for tool in minimap2_asm5 graphaligner vg_giraffe minigraph; do
            record_failure "$N" "$tool" "skipped" "linearise-and-map infeasible at n=$N"
        done
    fi
done

# =============================================================================
# PART 2 — HAPLOTYPE-PANEL SWEEP (Exp 3b)
# =============================================================================
log "=== PART 2: haplotype-panel sweep (fixed n=$HAP_N) ==="
HAPRESULTS="$OUTDIR/haplotype_scaling.tsv"
tsv_init "$HAPRESULTS" "n_snarls\thaplotypes\tsystem\trun\twall_sec\tpeak_rss_kb\tcpu_pct\texit"

for H in $HAP_PANELS; do
    log "----- panel h=$H (n=$HAP_N) -----"
    GA="$GRAPHDIR/ga_h${H}_n${HAP_N}.gfa"; GB="$GRAPHDIR/gb_h${H}_n${HAP_N}.gfa"
    if [[ ! -s "$GA" || ! -s "$GB" ]]; then
        python3 "$GEN" --mode haplotypes --n_snarls "$HAP_N" --haplotypes "$H" \
            --out_dir "$GRAPHDIR" --seed 42
    fi
    echo "$(sha_of "$GA")  $(basename "$GA")" >> "$MANIFEST"

    for RUN in $(seq 1 "$RUNS"); do
        m=$(measure_full "$G2G" "$GA" "$GB" --gaf /dev/null --threads "$NCPU")
        read -r wall peak user sys cpu ec <<< "$m"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$HAP_N" "$H" "g2g_parallel_${NCPU}t" "$RUN" "$wall" "$peak" "$cpu" "$ec" >> "$HAPRESULTS"
        log "  [g2g] h=$H run=$RUN wall=${wall}s rss=${peak}kB"
    done

    if [[ $HAVE_MINIMAP2 -eq 1 ]]; then
        FA="$GRAPHDIR/seqs_a_h${H}.fa"; FB="$GRAPHDIR/seqs_b_h${H}.fa"
        if timeout "$BASELINE_TIMEOUT" bash -c "$(declare -f linearise); linearise '$GA' '$FA'" 2>/dev/null; then
            linearise "$GB" "$FB" 2>/dev/null || true
            m=$(measure_full timeout "$BASELINE_TIMEOUT" \
                minimap2 -x asm5 -t "$NCPU" --secondary=no "$FA" "$FB" -o /dev/null)
            read -r wall peak user sys cpu ec <<< "$m"
            if [[ "$ec" == "0" ]]; then
                printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                    "$HAP_N" "$H" "minimap2_asm5" "1" "$wall" "$peak" "$cpu" "$ec" >> "$HAPRESULTS"
                log "  [minimap2] h=$H wall=${wall}s rss=${peak}kB"
            else
                record_failure "$HAP_N" "minimap2_asm5_h${H}" "timeout_or_error" "exit=$ec"
            fi
            rm -f "$FA" "$FB"
        fi
    fi
done

# =============================================================================
# PART 3 — SUMMARIES (mean +/- stderr)
# =============================================================================
log "=== PART 3: computing summaries ==="
summarise() {
    local in="$1" out="$2" keycols="$3"
    python3 - "$in" "$out" "$keycols" <<'PYEOF'
import sys, csv, math
from collections import defaultdict
infile, outfile, keycols = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
rows = defaultdict(lambda: {"wall": [], "rss": []})
with open(infile) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        if row.get("exit", "0") not in ("0", ""):
            continue
        key = tuple(row[k] for k in keycols)
        try:
            rows[key]["wall"].append(float(row["wall_sec"]))
            rows[key]["rss"].append(float(row["peak_rss_kb"]))
        except (ValueError, KeyError):
            pass
def stats(v):
    if not v: return (0.0, 0.0)
    m = sum(v)/len(v)
    if len(v) > 1:
        var = sum((x-m)**2 for x in v)/(len(v)-1)
        return (m, math.sqrt(var/len(v)))
    return (m, 0.0)
with open(outfile, "w") as f:
    hdr = keycols + ["mean_wall_sec","stderr_wall_sec","mean_peak_rss_gb","n_runs"]
    f.write("\t".join(hdr)+"\n")
    for key in sorted(rows):
        wm, ws = stats(rows[key]["wall"])
        rm, _  = stats(rows[key]["rss"])
        f.write("\t".join(list(key)+[f"{wm:.4f}",f"{ws:.4f}",f"{rm/1024/1024:.3f}",
                                     str(len(rows[key]['wall']))])+"\n")
print(f"[summary] wrote {outfile}")
PYEOF
}
summarise "$RESULTS"    "$OUTDIR/throughput_scaling_summary.tsv" "n_snarls,system"
summarise "$HAPRESULTS" "$OUTDIR/haplotype_scaling_summary.tsv"  "haplotypes,system"

log "=== Exp3 DONE ==="
log "  throughput:  $RESULTS"
log "  summary:     $OUTDIR/throughput_scaling_summary.tsv"
log "  failures:    $FAILMODE"
log "  haplotypes:  $HAPRESULTS"
log "  manifest:    $MANIFEST"
