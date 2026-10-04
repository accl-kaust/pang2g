#!/usr/bin/env bash
# ============================================================
# Experiment 6 — download ALL datasets on the LOGIN NODE.
#
# The compute nodes have NO internet access, so the
# Exp6 downloads cannot run as a Slurm job — they must run on a login node,
# which does have outbound network. Because the HPRC pull is large
# (~400 GB) and long-running, run this inside a tmux session so it survives
# disconnects:
#
#   tmux new -s g2g_dl
#   conda activate g2g_eval
#   bash experiments/download/download_all_login.sh
#   # detach: Ctrl-b then d        reattach: tmux attach -t g2g_dl
#
# This fetches, into $G2G_DATA_ROOT (default $HOME/g2g_data):
#   1. HPRC assemblies (AFR n=20 + EAS n=20, diploid)   ~400 GB
#   2. gnomAD-SV v4.1 sites VCF (+ index)
#   3. ClinVar GRCh38 VCF (+ index) and hg38.chrom.sizes
#
# Once this finishes, the compute-only stages run as Slurm jobs:
#   sbatch experiments/exp6_hprc/slurm_build_pangenomes.sh
#   sbatch experiments/exp6_hprc/slurm.sh
#
# Environment variables:
#   G2G_DATA_ROOT  base directory for data   [default: $HOME/g2g_data]
#   N_PER_POP      HPRC samples per pop      [default: 20]
# ============================================================
set -euo pipefail

G2G_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export G2G_DATA_ROOT="${G2G_DATA_ROOT:-$HOME/g2g_data}"
export N_PER_POP="${N_PER_POP:-20}"

# Run this on a LOGIN node: compute nodes usually have no outbound network.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    echo "[dl] WARNING: running inside a Slurm job; compute nodes have no internet." >&2
    echo "[dl] Run this on a login node inside tmux. Continuing in 5s..." >&2
    sleep 5
fi

# Preflight: the download/extraction toolchain must be on PATH. The most
# common failure is running in the (base) conda env instead of g2g_eval, which
# makes bcftools/tabix vanish and silently yields empty BED files.
MISSING=()
for tool in bcftools bgzip tabix wget python3; do
    command -v "$tool" &>/dev/null || MISSING+=("$tool")
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "[dl] ERROR: missing tools on PATH: ${MISSING[*]}" >&2
    echo "[dl]        Run 'conda activate g2g_eval' first." >&2
    exit 1
fi

echo "[dl] ============================================================"
echo "[dl] Exp6 dataset download (login node)"
echo "[dl] Host:       $HOST"
echo "[dl] Data root:  $G2G_DATA_ROOT"
echo "[dl] Started:    $(date)"
echo "[dl] ============================================================"

run_stage() {
    local name="$1"; shift
    echo ""
    echo "[dl] ── $name ─────────────────────────────────────────────"
    local start; start=$(date +%s)
    if bash "$@"; then
        echo "[dl] $name OK ($(( $(date +%s) - start ))s)"
    else
        echo "[dl] ERROR: $name failed (rc=$?). Fix and re-run; completed" >&2
        echo "[dl]        stages are skipped on re-run if data already exists." >&2
        exit 1
    fi
}

run_stage "HPRC assemblies" "$G2G_ROOT/experiments/download/download_hprc_assemblies.sh"
run_stage "GRCh38 reference" "$G2G_ROOT/experiments/download/download_grch38.sh"
run_stage "gnomAD-SV v4.1"  "$G2G_ROOT/experiments/download/download_gnomad_sv.sh"
run_stage "ClinVar (GRCh38)" "$G2G_ROOT/experiments/download/download_clinvar.sh"

echo ""
echo "[dl] ============================================================"
echo "[dl] All Exp6 downloads complete at $(date)"
echo "[dl] Next (compute nodes, as Slurm jobs):"
echo "[dl]   sbatch experiments/exp6_hprc/slurm_build_pangenomes.sh"
echo "[dl]   sbatch experiments/exp6_hprc/slurm.sh"
echo "[dl] ============================================================"
