#!/bin/bash
#SBATCH --job-name=bin_dastool
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs_binning/dastool/bin_dastool_%A_%a.out

# ============================================================
# 4d_refine_dastool.sh
# Author: Nickole Villabona
# Bin refinement and dereplication with DAS Tool
# ============================================================
# Integrates bins from MetaBAT2, MaxBin2, and CONCOCT using
# DAS Tool to select the best set of non-redundant bins.
# Evaluates bin quality using single-copy genes and outputs
# a refined set of MAGs.
#
# Input:  /work_meta/05_binning/<sample>/metabat/*.fa
#         /work_meta/05_binning/<sample>/maxbin/*.fasta
#         /work_meta/05_binning/<sample>/concoct/bins/*.fa
#         /work_meta/03_assembly/<sample>/final.contigs.fa
# Output: /work_meta/05_binning/<sample>/DASTool/
#         - <sample>_DASTool_DASTool_bins/*.fa (refined bins)
#         - <sample>_DASTool_summary.txt
#
# Tools: DAS Tool (diamond search engine, writes final bins)
#
# Usage: sbatch --array=1-113 4d_refine_dastool.sh
# ============================================================

set -euo pipefail

echo "[$(date)] ==== Starting DAS Tool SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} ===="

# ======== BASE PATHS ========
WORK_BASE=/dfs10/hammert-lab/nvillabo/Ch2/work_meta
ASSEMBLY_DIR=${WORK_BASE}/03_assembly
BIN_DIR=${WORK_BASE}/05_binning
LOG_DIR=${WORK_BASE}/logs_binning/dastool

mkdir -p "${BIN_DIR}" "${LOG_DIR}"

# ======== ACTIVATE ENVIRONMENT WITH DAS TOOL ========
module load miniconda3/23.5.2
module load samtools/1.15.1

# Initialize conda in non-interactive job
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate binning_env

DASTOOL_BIN=/data/homezvol3/nvillabo/.conda/envs/binning_env/bin
FASTA2BIN=${DASTOOL_BIN}/Fasta_to_Contig2Bin.sh
DASTOOL=${DASTOOL_BIN}/DAS_Tool

# ======== GET SAMPLE FROM ARRAY INDEX ========
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${WORK_BASE}/preprocessing/samples.txt")

if [ -z "${SAMPLE}" ]; then
    echo "[$(date)] No sample for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}, exiting."
    exit 0
fi

echo "[$(date)] Starting DAS Tool for sample ${SAMPLE}"

CONTIGS=${ASSEMBLY_DIR}/${SAMPLE}/final.contigs.fa

if [ ! -f "${CONTIGS}" ]; then
    echo "[$(date)] No assembly for ${SAMPLE}: ${CONTIGS}. Skipping."
    exit 0
fi

BIN_SAMPLE_DIR=${BIN_DIR}/${SAMPLE}
DASTOOL_DIR=${BIN_SAMPLE_DIR}/DASTool
mkdir -p "${DASTOOL_DIR}"
cd "${DASTOOL_DIR}"

# ======== 1) GENERATE CONTIG2BIN FILES FOR EACH BINNER ========

INPUTS=()
LABELS=()

# --- MetaBAT2 ---
METABAT_DIR=${BIN_SAMPLE_DIR}/metabat
METABAT_TSV=${DASTOOL_DIR}/${SAMPLE}_metabat2_scaffolds2bin.tsv

if ls "${METABAT_DIR}"/*metabat*.fa 1> /dev/null 2>&1; then
    echo "[$(date)] Generating contig2bin for MetaBAT2..."
    "${FASTA2BIN}" -e fa -i "${METABAT_DIR}" > "${METABAT_TSV}"
    if [ -s "${METABAT_TSV}" ]; then
        INPUTS+=("${METABAT_TSV}")
        LABELS+=("metabat2")
    else
        echo "[$(date)] WARNING: MetaBAT2 TSV is empty for ${SAMPLE}, will not use."
    fi
else
    echo "[$(date)] No MetaBAT2 bins found for ${SAMPLE}, skipping this binner."
fi

# --- MaxBin2 ---
MAXBIN_DIR=${BIN_SAMPLE_DIR}/maxbin
MAXBIN_TSV=${DASTOOL_DIR}/${SAMPLE}_maxbin2_scaffolds2bin.tsv

if ls "${MAXBIN_DIR}"/*.fasta 1> /dev/null 2>&1; then
    echo "[$(date)] Generating contig2bin for MaxBin2..."
    "${FASTA2BIN}" -e fasta -i "${MAXBIN_DIR}" > "${MAXBIN_TSV}"
    if [ -s "${MAXBIN_TSV}" ]; then
        INPUTS+=("${MAXBIN_TSV}")
        LABELS+=("maxbin2")
    else
        echo "[$(date)] WARNING: MaxBin2 TSV is empty for ${SAMPLE}, will not use."
    fi
else
    echo "[$(date)] No MaxBin2 bins found for ${SAMPLE}, skipping this binner."
fi

# --- CONCOCT ---
CONCOCT_DIR=${BIN_SAMPLE_DIR}/concoct/bins
CONCOCT_TSV=${DASTOOL_DIR}/${SAMPLE}_concoct_scaffolds2bin.tsv

if ls "${CONCOCT_DIR}"/*.fa 1> /dev/null 2>&1; then
    echo "[$(date)] Generating contig2bin for CONCOCT..."
    "${FASTA2BIN}" -e fa -i "${CONCOCT_DIR}" > "${CONCOCT_TSV}"
    if [ -s "${CONCOCT_TSV}" ]; then
        INPUTS+=("${CONCOCT_TSV}")
        LABELS+=("concoct")
    else
        echo "[$(date)] WARNING: CONCOCT TSV is empty for ${SAMPLE}, will not use."
    fi
else
    echo "[$(date)] No CONCOCT bins found for ${SAMPLE}, skipping this binner."
fi

# ======== 2) VERIFY WE HAVE AT LEAST ONE BINNER ========
if [ "${#INPUTS[@]}" -lt 1 ]; then
    echo "[$(date)] ERROR: No binners produced usable bins for ${SAMPLE}. Cannot run DAS Tool. Exiting."
    exit 0
fi

# Build comma-separated strings for DAS Tool
INPUT_STR=$(IFS=,; echo "${INPUTS[*]}")
LABEL_STR=$(IFS=,; echo "${LABELS[*]}")

echo "[$(date)] Using binners: ${LABEL_STR}"
echo "[$(date)] Contig2bin files: ${INPUT_STR}"

# ======== 3) RUN DAS TOOL ========

DASTOOL_PREFIX=${DASTOOL_DIR}/${SAMPLE}_DASTool

if [ -f "${DASTOOL_PREFIX}.DASTool_summary.txt" ]; then
    echo "[$(date)] DAS Tool already ran for ${SAMPLE}, skipping."
    exit 0
fi

echo "[$(date)] Running DAS Tool..."
"${DASTOOL}" \
  -i "${INPUT_STR}" \
  -l "${LABEL_STR}" \
  -c "${CONTIGS}" \
  -o "${DASTOOL_PREFIX}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --search_engine diamond \
  --write_bins

echo "[$(date)] Finished DAS Tool for ${SAMPLE}"
