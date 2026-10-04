#!/bin/bash --login
# Camera-ready re-measurement, part 3 (real chr22 data already staged under the
# previous checkout's demo/data): throughput vs. linearize-and-map (demo/10)
# and the Q4 AFR-vs-EAS run (demo/09), both with the superbubble decomposer.
#     DATA=~/Research/G2G/G2G_Project_Summer_2026/demo/data sbatch experiments/accl/cr_realdata.slurm
#SBATCH --job-name=g2g_cr_real
#SBATCH --partition=build
#SBATCH --cpus-per-task=48
#SBATCH --mem=192G
#SBATCH --time=10:00:00
#SBATCH --output=logs/cr_real_%j.out
#SBATCH --error=logs/cr_real_%j.err
set -uo pipefail
R="${SLURM_SUBMIT_DIR:-$PWD}"; cd "$R"; mkdir -p logs results/accl
DATA="${DATA:-$HOME/Research/G2G/G2G_Project_Summer_2026/demo/data}"
OLD="$HOME/Research/G2G/G2G_Project_Summer_2026"
# python3 (+numpy) from conda env g2g; modern htslib tools from env g2g_bio; vg and
# minimap2 from the previous checkout's demo/tools.  The demo scripts honour these
# overrides (VG=, BCFTOOLS=, TABIX=, SAMTOOLS=, MINIMAP2=).
export PATH=$HOME/miniconda3/envs/g2g/bin:$HOME/miniconda3/bin:/opt/slurm/current/bin:$PATH
BIO=$HOME/miniconda3/envs/g2g_bio/bin
export VG=$OLD/demo/tools/vg BCFTOOLS=$BIO/bcftools TABIX=$BIO/tabix SAMTOOLS=$BIO/samtools MINIMAP2=$OLD/demo/tools/minimap2
echo "[cr] $(date) host=$(hostname) vg=$($VG version 2>&1 | head -1) bcftools=$($BCFTOOLS --version | head -1) minimap2=$($MINIMAP2 --version) python=$(python3 --version)"
for i in $(seq 1 240); do [ -x build/release/g2g ] && break; sleep 30; done
[ -x build/release/g2g ] || { echo "no binary"; exit 1; }
ln -sfn release/g2g build/g2g
python3 -c "import numpy" 2>/dev/null || { echo "[cr] numpy missing in PATH python; enrichment will use the slow stdlib path"; }

echo "[cr] === Q4: AFR vs EAS chr22 (demo/09 real mode) ==="
REF_FA="$DATA/human_g1k_v37.fasta" PANEL_TMPL="$DATA/ALL.chr{}.phase3.genotypes.vcf.gz" \
SAMPLES_A="$DATA/afr.samples" SAMPLES_B="$DATA/eas.samples" \
GNOMAD_SV_VCF="$DATA/gnomad_v2.1_sv.sites.vcf.gz" CLINVAR_VCF="$DATA/clinvar_GRCh37.vcf.gz" \
CHROM_SIZES="$DATA/b37.chrom.sizes" CHROMS="22" AFR_FIELD=AFR_AF EAS_FIELD=EAS_AF \
bash demo/09_hprc_biology.sh 2>&1 | grep -v "^snarl" | tee logs/09_real_chr22_accl.log | tail -12
cp demo/work/biology/biology_results.tsv results/accl/biology_results.tsv 2>/dev/null
python3 demo/hprc_enrich.py --delta demo/work/biology/real/delta.vcf --gnomad-sv "$DATA/gnomad_v2.1_sv.sites.vcf.gz" \
    --sizes "$DATA/b37.chrom.sizes" --afr-field AFR_AF --eas-field EAS_AF --window 400 --perms 1000 \
    --fst-sweep 0.05,0.1,0.2,0.3,0.5 --out results/accl/fst_sweep_chr22.tsv 2>&1 | tail -7

echo "[cr] === throughput vs linearize-and-map (demo/10, h sweep) ==="
REF_FA="$DATA/human_g1k_v37.fasta" PANEL_TMPL="$DATA/ALL.chr{}.phase3.genotypes.vcf.gz" \
SAMPLES_A="$DATA/afr.samples" SAMPLES_B="$DATA/eas.samples" CONTIG=22 H_SWEEP="2 5 10 25 50 100 200" \
bash demo/10_throughput_baseline.sh 2>&1 | tee logs/10_throughput_accl.log | tail -14
cp demo/work/throughput/throughput.tsv results/accl/throughput_chr22.tsv 2>/dev/null
echo "[cr] $(date) DONE_REAL"
