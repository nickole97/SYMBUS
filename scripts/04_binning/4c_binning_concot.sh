#!/bin/bash
#SBATCH --job-name=bin_concoct
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs_binning/concoct/bin_concoct_%A_%a.out

# ============================================================
# 4c_binning_concoct.sh
# Author: Nickole Villabona
# Genome binning with CONCOCT
# ============================================================
# Bins assembled contigs into MAGs using CONCOCT, which uses
# sequence composition and coverage across multiple samples.
# Contigs are first cut into 10kb chunks before clustering.
#
# Input:  /work_meta/03_assembly/<sample>/final.contigs.fa
#         /work_meta/04_mapping/<sample>/<sample>.sorted.bam
# Output: /work_meta/05_binning/<sample>/concoct/bins/
#         - *.fa (one file per bin)
#
# Tools: CONCOCT (cut contigs at 10kb, merge after clustering)
#
# Usage: sbatch --array=1-113%7 4c_binning_concoct.sh
# ============================================================

set -euo pipefail

echo "[$(date)] ==== Starting CONCOCT SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} ===="

# ======== BASE PATHS ========
WORK_BASE=/dfs10/hammert-lab/nvillabo/Ch2/work_meta
ASSEMBLY_DIR=${WORK_BASE}/03_assembly
MAP_DIR=${WORK_BASE}/04_mapping
BIN_DIR=${WORK_BASE}/05_binning
LOG_DIR=${WORK_BASE}/logs_binning/concoct

mkdir -p "${BIN_DIR}" "${LOG_DIR}"

# ======== LOAD SAMTOOLS FOR CONCOCT ========
module load samtools/1.15.1

# ======== GET SAMPLE FROM ARRAY INDEX ========
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${WORK_BASE}/preprocessing/samples.txt")

if [ -z "${SAMPLE}" ]; then
    echo "[$(date)] No sample for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}, exiting."
    exit 0
fi

echo "[$(date)] Starting CONCOCT for sample ${SAMPLE}"

CONTIGS=${ASSEMBLY_DIR}/${SAMPLE}/final.contigs.fa
BAM=${MAP_DIR}/${SAMPLE}/${SAMPLE}.sorted.bam

# ======== BASIC CHECKS ========
if [ ! -f "${CONTIGS}" ]; then
    echo "[$(date)] No assembly for ${SAMPLE}: ${CONTIGS}. Skipping."
    exit 0
fi

if [ ! -f "${BAM}" ]; then
    echo "[$(date)] No BAM file for ${SAMPLE}: ${BAM}. Skipping."
    exit 0
fi

# ======== PATHS TO CONCOCT TOOLS ========
CONCOCT_ENV=/data/homezvol3/nvillabo/.conda/envs/binning_env/bin
CUTUP=${CONCOCT_ENV}/cut_up_fasta.py
COVTAB=${CONCOCT_ENV}/concoct_coverage_table.py
CONCOCT=${CONCOCT_ENV}/concoct
MERGE=${CONCOCT_ENV}/merge_cutup_clustering.py
EXTRACT=${CONCOCT_ENV}/extract_fasta_bins.py

BIN_SAMPLE_DIR=${BIN_DIR}/${SAMPLE}/concoct
mkdir -p "${BIN_SAMPLE_DIR}"
cd "${BIN_SAMPLE_DIR}"

echo "[$(date)] Working in ${BIN_SAMPLE_DIR}"

# ======== 1) CUT CONTIGS INTO 10KB CHUNKS ========
CUT_FASTA=contigs_10K.fa
CUT_BED=contigs_10K.bed

if [ ! -s "${CUT_FASTA}" ] || [ ! -s "${CUT_BED}" ]; then
    echo "[$(date)] Cutting contigs into 10kb chunks..."
    "${CUTUP}" "${CONTIGS}" -c 10000 -o 0 --merge_last -b "${CUT_BED}" > "${CUT_FASTA}"
else
    echo "[$(date)] Cut files already exist, skipping cut_up_fasta."
fi

# ======== 2) COVERAGE TABLE ========
COV_TABLE=coverage_table.tsv

if [ ! -s "${COV_TABLE}" ]; then
    echo "[$(date)] Calculating coverage table..."
    "${COVTAB}" "${CUT_BED}" "${BAM}" > "${COV_TABLE}"
else
    echo "[$(date)] Coverage table already exists, skipping."
fi

# ======== 3) RUN CONCOCT ON CUT CONTIGS ========
CONCOCT_PREFIX=concoct_out
CLUSTER_RAW=${CONCOCT_PREFIX}_clustering_gt1000.csv

if [ ! -s "${CLUSTER_RAW}" ]; then
    echo "[$(date)] Running CONCOCT..."
    "${CONCOCT}" \
      --composition_file "${CUT_FASTA}" \
      --coverage_file "${COV_TABLE}" \
      --threads "${SLURM_CPUS_PER_TASK}" \
      -b "${CONCOCT_PREFIX}"
else
    echo "[$(date)] Raw CONCOCT clustering already exists (${CLUSTER_RAW}), skipping."
fi

if [ ! -s "${CLUSTER_RAW}" ]; then
    echo "[$(date)] ERROR: ${CLUSTER_RAW} not generated. Cannot continue. Exiting."
    exit 1
fi

# ======== 4) MERGE CLUSTERING FROM CHUNKS TO ORIGINAL CONTIGS ========
MERGED_CLUSTERING=${CONCOCT_PREFIX}_clustering_merged.csv

if [ ! -s "${MERGED_CLUSTERING}" ]; then
    echo "[$(date)] Merging clustering from chunks to original contigs..."
    "${MERGE}" "${CLUSTER_RAW}" > "${MERGED_CLUSTERING}"
else
    echo "[$(date)] Merged clustering already exists (${MERGED_CLUSTERING}), skipping merge."
fi

if [ ! -s "${MERGED_CLUSTERING}" ]; then
    echo "[$(date)] ERROR: ${MERGED_CLUSTERING} not generated. Cannot extract bins. Exiting."
    exit 1
fi

# ======== 5) EXTRACT BINS TO FASTA ========
BINS_DIR=${BIN_SAMPLE_DIR}/bins
mkdir -p "${BINS_DIR}"

# Skip if bins already exist
if ls "${BINS_DIR}"/*.fa 1> /dev/null 2>&1; then
    echo "[$(date)] CONCOCT bins already exist for ${SAMPLE} in ${BINS_DIR}, exiting."
    exit 0
fi

echo "[$(date)] Extracting CONCOCT bins..."
"${EXTRACT}" \
  "${CONTIGS}" \
  "${MERGED_CLUSTERING}" \
  --output_path "${BINS_DIR}"

echo "[$(date)] Finished CONCOCT for ${SAMPLE}"
