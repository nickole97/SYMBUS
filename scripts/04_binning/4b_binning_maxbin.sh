#!/bin/bash
#SBATCH --job-name=bin_maxbin
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs_binning/bin_maxbin_%A_%a.out

# ============================================================
# 4b_binning_maxbin.sh
# Author: Nickole Villabona
# Genome binning with MaxBin2
# ============================================================
# Bins assembled contigs into MAGs using MaxBin2, which uses
# tetranucleotide frequency and abundance information from
# read mapping depth.
#
# Input:  /work_meta/03_assembly/<sample>/final.contigs.fa
#         /work_meta/04_mapping/<sample>/<sample>_depth.txt
# Output: /work_meta/05_binning/<sample>/maxbin/
#         - <sample>_maxbin.*.fasta (one file per bin)
#
# Tools: MaxBin2 (min-contig-len=1500)
#
# Usage: sbatch --array=1-113%5 4b_binning_maxbin.sh
# ============================================================

set -euo pipefail

WORK_BASE=/dfs10/hammert-lab/nvillabo/Ch2/work_meta
ASSEMBLY_DIR=${WORK_BASE}/03_assembly
MAP_DIR=${WORK_BASE}/04_mapping
BIN_DIR=${WORK_BASE}/05_binning

mkdir -p "${BIN_DIR}" "${WORK_BASE}/logs_binning/maxbin"

# ======== GET SAMPLE FROM ARRAY INDEX ========
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${WORK_BASE}/preprocessing/samples.txt")

if [ -z "${SAMPLE}" ]; then
    echo "[$(date)] No sample for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}, exiting."
    exit 0
fi

echo "[$(date)] Starting MaxBin2 for sample ${SAMPLE}"

CONTIGS=${ASSEMBLY_DIR}/${SAMPLE}/final.contigs.fa
DEPTH=${MAP_DIR}/${SAMPLE}/${SAMPLE}_depth.txt

# ======== BASIC CHECKS ========
if [ ! -f "${CONTIGS}" ]; then
    echo "[$(date)] No assembly for ${SAMPLE}: ${CONTIGS}. Skipping."
    exit 0
fi

if [ ! -f "${DEPTH}" ]; then
    echo "[$(date)] No depth file for ${SAMPLE}: ${DEPTH}. Skipping."
    exit 0
fi

# ======== PATH TO MAXBIN2 ========
MAXBIN=~/.conda/envs/binning_env/bin/run_MaxBin.pl

BIN_SAMPLE_DIR=${BIN_DIR}/${SAMPLE}/maxbin
mkdir -p "${BIN_SAMPLE_DIR}"

# Skip if bins already exist
if ls "${BIN_SAMPLE_DIR}/${SAMPLE}_maxbin."*.fasta 1> /dev/null 2>&1; then
    echo "[$(date)] MaxBin2 bins already exist for ${SAMPLE} in ${BIN_SAMPLE_DIR}, exiting."
    exit 0
fi

cd "${BIN_SAMPLE_DIR}"

echo "[$(date)] Running MaxBin2..."
"${MAXBIN}" \
  -contig "${CONTIGS}" \
  -abund "${DEPTH}" \
  -out "${SAMPLE}_maxbin" \
  -thread "${SLURM_CPUS_PER_TASK}" \
  -min_contig_length 1500

echo "[$(date)] Finished MaxBin2 for ${SAMPLE}"
