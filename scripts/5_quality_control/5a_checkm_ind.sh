#!/bin/bash
#SBATCH --job-name=checkm_bins
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs_binning/checkm/checkm_%A_%a.out

# ============================================================
# 5a_checkm_individual.sh
# Author: Nickole Villabona
# MAG quality assessment with CheckM
# ============================================================
# Assesses the quality of bins from DAS Tool using CheckM.
# Calculates completeness, contamination, and strain heterogeneity
# for each MAG using lineage-specific marker genes.
#
# Input:  /work_meta/05_binning/<sample>/DASTool/<sample>_DASTool_DASTool_bins/*.fa
# Output: /work_meta/05_binning/<sample>/CheckM/
#         - storage/bin_stats_ext.tsv (detailed stats)
#         - lineage.ms (marker set used)
#
# Tools: CheckM lineage_wf (lineage-specific workflow)
#
# Usage: sbatch --array=1-113 5a_checkm_individual.sh
# Runtime: ~6-12 hours per sample depending on number of bins
# ============================================================

set -euo pipefail

echo "[$(date)] ==== Starting CheckM SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} ===="

# ======== BASE PATHS ========
WORK_BASE=/dfs10/hammert-lab/nvillabo/Ch2/work_meta
BIN_DIR=${WORK_BASE}/05_binning
LOG_DIR=${WORK_BASE}/logs_binning/checkm

mkdir -p "${LOG_DIR}"

# ======== DEFINE CHECKM AND DATABASE ========
module load miniconda3/23.5.2

CHECKM_BIN=/data/homezvol3/nvillabo/.conda/envs/checkm_env/bin/checkm
export CHECKM_DATA_PATH=/data/homezvol3/nvillabo/.conda/envs/checkm_env/checkm_data

echo "[$(date)] Using CHECKM_BIN=${CHECKM_BIN}"
echo "[$(date)] CHECKM_DATA_PATH=${CHECKM_DATA_PATH}"

# ======== GET SAMPLE FROM ARRAY INDEX ========
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${WORK_BASE}/preprocessing/samples.txt")

if [ -z "${SAMPLE}" ]; then
    echo "[$(date)] No sample for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}, exiting."
    exit 0
fi

echo "[$(date)] Starting CheckM for sample ${SAMPLE}"

BIN_SAMPLE_DIR=${BIN_DIR}/${SAMPLE}
DASTOOL_BINS=${BIN_SAMPLE_DIR}/DASTool/${SAMPLE}_DASTool_DASTool_bins

if [ ! -d "${DASTOOL_BINS}" ]; then
    echo "[$(date)] No DAS Tool bins folder for ${SAMPLE}: ${DASTOOL_BINS}. Skipping."
    exit 0
fi

# Check if there are any MAGs (.fa files)
if ! ls "${DASTOOL_BINS}"/*.fa 1> /dev/null 2>&1; then
    echo "[$(date)] No MAGs (.fa) in ${DASTOOL_BINS}. Skipping."
    exit 0
fi

CHECKM_OUT=${BIN_SAMPLE_DIR}/CheckM
mkdir -p "${CHECKM_OUT}"

# Skip if already ran
if [ -s "${CHECKM_OUT}/storage/bin_stats_ext.tsv" ]; then
    echo "[$(date)] CheckM already ran for ${SAMPLE} (bin_stats_ext.tsv exists). Exiting."
    exit 0
fi

echo "[$(date)] Running checkm lineage_wf on ${DASTOOL_BINS}"
"${CHECKM_BIN}" lineage_wf \
  -t "${SLURM_CPUS_PER_TASK}" \
  -x fa \
  "${DASTOOL_BINS}" \
  "${CHECKM_OUT}"

echo "[$(date)] Finished CheckM for ${SAMPLE}"
