#!/bin/bash
#SBATCH --job-name=covMAG
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs/coverm/covMAG_%A_%a.out

# ============================================================
# 7a_coverm_genome.sh
# Author: Nickole Villabona
# Quantify MAG abundance across samples with CoverM
# ============================================================
# Maps non-host reads to dereplicated MAGs and calculates
# abundance metrics including mean coverage, covered fraction,
# and TPM (transcripts per million) for each MAG in each sample.
#
# NOTE: This approach may be replaced with inStrain following
#       Meline's suggestion for strain-level profiling.
#
# Input:  /work_meta/preprocessing/02_host_removal/<sample>/*_nonhost.fastq.gz
#         /work_meta/06_MAGs_selected/*.fa
# Output: /work_meta/07_abundance/<sample>_coverm_genome.tsv
#
# Tools: CoverM genome mode (minimap2-sr mapper)
#        Thresholds: 95% identity, 75% alignment
#
# Usage: sbatch --array=1-113 7a_coverm_genome.sh
# ============================================================

set -euo pipefail

echo "[$(date)] ==== Starting CoverM SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} ===="

# ======== BASE PATHS ========
WORK_BASE=/dfs10/hammert-lab/nvillabo/Ch2/work_meta
MAG_DIR=${WORK_BASE}/06_MAGs_selected
READ_DIR=${WORK_BASE}/preprocessing/02_host_removal
OUT_DIR=${WORK_BASE}/07_abundance
LOG_DIR=${WORK_BASE}/logs/coverm

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

# ======== ABSOLUTE PATHS TO TOOLS ========
COVERM=~/.conda/envs/coverm_env/bin/coverm
MINIMAP2=~/.conda/envs/coverm_env/bin/minimap2

# ======== GET SAMPLE FROM ARRAY INDEX ========
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${WORK_BASE}/preprocessing/samples.txt")

if [ -z "${SAMPLE}" ]; then
    echo "[$(date)] No sample for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}, exiting."
    exit 0
fi

echo "[$(date)] Processing sample ${SAMPLE}"

R1=${READ_DIR}/${SAMPLE}/${SAMPLE}_R1_nonhost.fastq.gz
R2=${READ_DIR}/${SAMPLE}/${SAMPLE}_R2_nonhost.fastq.gz

if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
    echo "[$(date)] Missing reads for ${SAMPLE}: ${R1} or ${R2}. Skipping."
    exit 0
fi

OUT_TSV=${OUT_DIR}/${SAMPLE}_coverm_genome.tsv

if [ -s "${OUT_TSV}" ]; then
    echo "[$(date)] ${OUT_TSV} already exists — skipping."
    exit 0
fi

echo "[$(date)] Running coverm genome..."

${COVERM} genome \
  -1 "${R1}" \
  -2 "${R2}" \
  --genome-fasta-directory "${MAG_DIR}" \
  --genome-fasta-extension fa \
  --mapper minimap2-sr \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --methods mean covered_fraction tpm \
  --min-read-percent-identity 95 \
  --min-read-aligned-percent 75 \
  -o "${OUT_TSV}"

echo "[$(date)] CoverM finished for ${SAMPLE}"
