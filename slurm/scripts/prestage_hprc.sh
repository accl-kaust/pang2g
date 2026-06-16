#!/bin/bash
# =============================================================================
# prestage_hprc.sh — Download HPRC assemblies + truth sets on the LOGIN NODE.
# =============================================================================
# WHY THIS EXISTS
#   The ROCS testbed compute nodes (mcnode*) sit on a private network reached
#   only via VPN; they have no verified outbound route to AWS S3. exp6's S1
#   download therefore fails INSTANTLY on a compute node ("connection refused",
#   ~1 s per file — exactly the failure in exp6-*.err). The fix is to fetch all
#   external data ONCE, here, on the login node (mcmgt01, which does have
#   outbound access), into the shared network filesystem. exp6 then runs on a
#   compute node and only *verifies + uses* the pre-staged files (never downloads).
#
#   This script also fixes the three URL bugs in the old exp6 fetch_assembly():
#     1. base path had a spurious '/assemblies'  ->  must end at '/working/HPRC'
#     2. missing the 'year1_f1_assembly_v2_genbank/' subdirectory
#     3. filename used 'f1_assembly_v2'  ->  must be 'f1_assembly_v2_genbank'
#   Rather than hand-build URLs (the original mistake), we drive every download
#   from the OFFICIAL HPRC index file, which is the single source of truth and
#   also carries the per-haplotype SHA-256 for real verification.
#
# RUN (interactively, on the login node, NOT via sbatch):
#     bash slurm/scripts/prestage_hprc.sh
#
# It is idempotent and resumable: already-present, sha-verified files are skipped.
#
# CONFIG (env overrides):
#   G2G_REPO        repo root (default: two levels up)
#   G2G_DATA_ROOT   where to stage data (default: $REPO/data/hprc) — point at the
#                   shared network FS so compute nodes can read it
#   G2G_AFR_LIST / G2G_EAS_LIST   panel sample-ID lists
#   HPRC_INDEX_URL  override the index URL (default: HPP_Year1 genbank v2 index)
# =============================================================================
set -euo pipefail

if [[ -n "${G2G_REPO:-}" ]]; then REPO="$G2G_REPO"
else REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; fi
SCRIPT_DIR="$REPO/slurm/scripts"

DATA_ROOT="${G2G_DATA_ROOT:-$REPO/data/hprc}"
ASM_DIR="$DATA_ROOT/assemblies"
REF_DIR="$DATA_ROOT/reference"
INDEX_DIR="$DATA_ROOT/index"
mkdir -p "$ASM_DIR" "$REF_DIR" "$INDEX_DIR"

AFR_SAMPLES="${G2G_AFR_LIST:-$SCRIPT_DIR/panels/afr_samples.txt}"
EAS_SAMPLES="${G2G_EAS_LIST:-$SCRIPT_DIR/panels/eas_samples.txt}"

# Official HPRC Year-1 freeze-2 GenBank assembly index (sample -> S3/GCP URLs + sha256).
HPRC_INDEX_URL="${HPRC_INDEX_URL:-https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index}"
INDEX_FILE="$INDEX_DIR/Year1_assemblies_v2_genbank.index"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
sha_of() { sha256sum "$1" 2>/dev/null | awk '{print $1}' || echo NA; }

# S3 https endpoint for the (no-egress) HPRC public bucket.
s3_to_https() { sed 's%^s3://human-pangenomics/%https://s3-us-west-2.amazonaws.com/human-pangenomics/%'; }

# ── Lazy network access ───────────────────────────────────────────────────────
# The preflight + index fetch are DEFERRED until something actually needs to be
# downloaded. If every assembly is already present + SHA-verified, the script
# runs entirely offline — no network check, no index download, no wget needed.
S3_HOST="s3-us-west-2.amazonaws.com"
_NETWORK_CHECKED=0   # 0=not yet, 1=reachable, 2=unreachable

# Called the first time we need to download. Checks wget presence, S3 reach,
# and fetches the index if absent. Exits with guidance on failure.
require_network() {
    [[ "$_NETWORK_CHECKED" -ne 0 ]] && return $(( _NETWORK_CHECKED == 1 ? 0 : 1 ))

    command -v wget >/dev/null 2>&1 \
        || { log "FATAL: wget not found. 'conda install -c conda-forge wget' then re-run."; _NETWORK_CHECKED=2; return 1; }

    log "preflight: checking outbound access to $S3_HOST ..."
    local ok=0
    if wget --timeout=15 --tries=1 --spider -q "https://$S3_HOST/" 2>/dev/null; then
        ok=1
    elif command -v curl >/dev/null 2>&1 && curl -sS --max-time 15 -I "https://$S3_HOST/" >/dev/null 2>&1; then
        ok=1
    fi

    if [[ "$ok" -ne 1 ]]; then
        log "FATAL: cannot reach $S3_HOST from this machine ($(hostname))."
        log "       The ROCS testbed login/compute nodes are VPN-gated and may have"
        log "       no outbound internet. You have two options:"
        log "  (A) Download the assemblies on a machine that HAS internet (your laptop),"
        log "      then copy them into the shared FS, e.g.:"
        log "        rsync -avP ./assemblies/  <user>@mcmgt01:$ASM_DIR/"
        log "      exp6 only needs the .fa.gz files sitting in $ASM_DIR (it verifies them)."
        log "  (B) If a proxy is available on the testbed, export it and re-run:"
        log "        export https_proxy=http://<proxy-host>:<port> http_proxy=\$https_proxy"
        log "       (ask the ROCS admins for the proxy address if unsure)."
        if [[ "${G2G_SKIP_PREFLIGHT:-0}" -eq 1 ]]; then
            log "G2G_SKIP_PREFLIGHT=1 set — continuing despite failed preflight."
        else
            _NETWORK_CHECKED=2; return 1
        fi
    else
        log "preflight OK: $S3_HOST reachable."
    fi
    _NETWORK_CHECKED=1
    return 0
}

# Ensure the HPRC index is present (fetches it if needed + network available).
require_index() {
    if [[ -s "$INDEX_FILE" ]] && head -1 "$INDEX_FILE" | grep -q 'hap1_aws_fasta'; then
        return 0   # already have it
    fi
    # Need to download it.
    if ! require_network; then
        log "cannot fetch HPRC index — network unavailable"
        return 1
    fi
    log "downloading HPRC index -> $INDEX_FILE"
    wget --tries=3 --waitretry=5 --timeout=30 --show-progress -O "$INDEX_FILE" "$HPRC_INDEX_URL" \
        || { log "FATAL: could not fetch HPRC index from $HPRC_INDEX_URL"; return 1; }
    head -1 "$INDEX_FILE" | grep -q 'hap1_aws_fasta' \
        || { log "FATAL: $INDEX_FILE does not look like the HPRC index"; return 1; }
    return 0
}

# Try to load the index (needed for resolving filenames + SHAs).
if require_index; then
    log "index ready ($(wc -l < "$INDEX_FILE") rows)"
    INDEX_OK=1
else
    log "HPRC index not available — can only verify already-downloaded files by filename"
    INDEX_OK=0
fi

# ── Resolver: sample+hap -> (https_url, expected_sha) from the index ─────────
# paternal => hap1 (url col 2, sha col 6); maternal => hap2 (url col 3, sha col 7)
resolve() {
    local sample="$1" hap="$2" urlcol shacol
    case "$hap" in
        paternal) urlcol=2; shacol=6 ;;
        maternal) urlcol=3; shacol=7 ;;
        *) return 2 ;;
    esac
    awk -F'\t' -v s="$sample" -v uc="$urlcol" -v sc="$shacol" '
        $1==s { print $uc "\t" $sc; found=1 }
        END   { if(!found) exit 3 }' "$INDEX_FILE"
}

# ── Fetch/verify one haplotype assembly ───────────────────────────────────────
# Priority order:
#   1. If index available: resolve canonical filename + SHA, check if already
#      present + verified → skip.  If missing → download (triggers preflight).
#   2. If no index: look for any file matching the expected naming pattern in
#      $ASM_DIR → report found (unverified) or missing.
fetch_one() {
    local panel="$1" sample="$2" hap="$3"

    # ── With index: full resolve + SHA verify ─────────────────────────────────
    if [[ "$INDEX_OK" -eq 1 ]]; then
        local row s3url sha
        if ! row=$(resolve "$sample" "$hap"); then
            log "  [$panel] $sample.$hap NOT in index — check the sample ID"; return 1
        fi
        s3url=$(printf '%s' "$row" | cut -f1)
        sha=$(  printf '%s' "$row" | cut -f2)
        local url; url=$(printf '%s' "$s3url" | s3_to_https)
        local fname; fname=$(basename "$url")
        local dest="$ASM_DIR/$fname"

        # Already present + SHA verified → skip (no network needed).
        if [[ -s "$dest" ]] && [[ "$(sha_of "$dest")" == "$sha" ]]; then
            log "  [$panel] $fname already present + verified (skip)"
            echo "$dest"; return 0
        fi

        # Need to download — check network first (lazy, runs only once).
        if ! require_network; then
            if [[ -s "$dest" ]]; then
                log "  [$panel] $fname present but SHA unverified (network unavailable for re-download)"
                echo "$dest"; return 0   # use what we have
            fi
            log "  [$panel] $fname MISSING and cannot download (no network)"
            return 1
        fi

        log "  [$panel] fetching $fname"
        if ! wget --continue --tries=3 --waitretry=15 --timeout=60 --show-progress \
                  -O "$dest" "$url"; then
            log "  [$panel] download FAILED for $fname (network/route problem?)"
            rm -f "$dest"
            return 1
        fi

        local got; got=$(sha_of "$dest")
        if [[ "$got" != "$sha" ]]; then
            log "  [$panel] SHA MISMATCH for $fname (got $got want $sha) — deleting"
            rm -f "$dest"; return 1
        fi
        log "  [$panel] OK $fname (sha verified)"
        echo "$dest"; return 0
    fi

    # ── Without index: best-effort filename match ─────────────────────────────
    # Try the canonical HPRC naming pattern.
    local fname="${sample}.${hap}.f1_assembly_v2_genbank.fa.gz"
    local dest="$ASM_DIR/$fname"
    if [[ -s "$dest" ]]; then
        log "  [$panel] $fname found (no index — SHA not verified)"
        echo "$dest"; return 0
    fi
    log "  [$panel] $fname MISSING (and no index to resolve the download URL)"
    return 1
}

stage_panel() {
    local panel="$1" listfile="$2" ok=0 fail=0
    [[ -s "$listfile" ]] || { log "panel list $listfile missing/empty"; return 1; }
    while IFS= read -r sample; do
        [[ -z "$sample" || "$sample" == \#* ]] && continue
        for hap in paternal maternal; do
            if fetch_one "$panel" "$sample" "$hap" >/dev/null; then ok=$((ok+1)); else fail=$((fail+1)); fi
        done
    done < "$listfile"
    log "$panel: $ok ok, $fail failed"
    [[ "$fail" -eq 0 ]]
}

log "=== Pre-staging HPRC assemblies (login node) ==="
RC=0
stage_panel AFR "$AFR_SAMPLES" || RC=1
stage_panel EAS "$EAS_SAMPLES" || RC=1

# ── Truth sets reminder (cannot be auto-fetched without stable public URLs) ──
for f in "GRCh38 FASTA:$REF_DIR/GRCh38.fa" \
         "gnomAD-SV v4 BED:$REF_DIR/gnomad_sv_v4.bed" \
         "ClinVar SV BED:$REF_DIR/clinvar_sv.bed"; do
    name="${f%%:*}"; path="${f##*:}"
    if [[ -s "$path" ]]; then log "truth set present: $name -> $path"
    else log "MISSING truth set: $name — stage it at $path (see README)"; fi
done

log "=== Pre-stage complete (rc=$RC). Assemblies in: $ASM_DIR ==="
log "Now submit exp6 (it will verify + use these files, no download):"
log "  sbatch slurm/scripts/exp6_biological.sh"
exit "$RC"