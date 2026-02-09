#!/bin/bash
#SBATCH --job-name=bin_metabat
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs_binning/bin_metabat_%A_%a.out

# ============================================================
# 4a_binning_metabat.sh
# Author: Nickole Villabona
# Genome binning with MetaBAT2
# ============================================================
# Bins assembled contigs into MAGs (Metagenome-Assembled Genomes)
# using MetaBAT2, which uses tetranucleotide frequency and 
# coverage depth information.
#
# Input:  /work_meta/03_assembly/<sample>/final.contigs.fa
#         /work_meta/04_mapping/<sample>/<sample>_depth.txt
# Output: /work_meta/05_binning/<sample>/metabat/
#         - <sample>_metabat.*.fa (one file per bin)
#
# Tools: MetaBAT2 (min-contig-len=1500)
#
# Usage: sbatch --array=1-113%5 4a_binning_metabat.sh
# ============================================================

set -euo pipefail

WORK_BASE=/dfs10/hammert-lab/nvillabo/Ch2/work_meta
ASSEMBLY_DIR=${WORK_BASE}/03_assembly
MAP_DIR=${WORK_BASE}/04_mapping
BIN_DIR=${WORK_BASE}/05_binning

mkdir -p "${BIN_DIR}" "${WORK_BASE}/logs_binning"

# ======== GET SAMPLE FROM ARRAY INDEX ========
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${WORK_BASE}/preprocessing/samples.txt")

if [ -z "${SAMPLE}" ]; then
    echo "[$(date)] No sample for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}, exiting."
    exit 0
fi

echo "[$(date)] Starting MetaBAT2 for sample ${SAMPLE}"

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

# ======== PATH TO METABAT2 ========
METABAT=~/.conda/envs/binning_env/bin/metabat2

BIN_SAMPLE_DIR=${BIN_DIR}/${SAMPLE}/metabat
mkdir -p "${BIN_SAMPLE_DIR}"

# Skip if bins already exist
if ls "${BIN_SAMPLE_DIR}/${SAMPLE}_metabat."*.fa 1> /dev/null 2>&1; then
    echo "[$(date)] MetaBAT2 bins already exist for ${SAMPLE} in ${BIN_SAMPLE_DIR}, exiting."
    exit 0
fi

echo "[$(date)] Running metabat2..."
"${METABAT}" \
  -i "${CONTIGS}" \
  -a "${DEPTH}" \
  -o "${BIN_SAMPLE_DIR}/${SAMPLE}_metabat" \
  -m 1500 \
  -t "${SLURM_CPUS_PER_TASK}"

echo "[$(date)] Finished MetaBAT2 for ${SAMPLE}"
