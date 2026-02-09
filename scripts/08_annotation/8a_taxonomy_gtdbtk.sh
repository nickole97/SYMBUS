#!/bin/bash
#SBATCH --job-name=gtdbtk
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs/gtdbtk/gtdbtk_%A.out

# ============================================================
# 8a_taxonomy_gtdbtk.sh
# Author: Nickole Villabona
# Taxonomic classification with GTDB-Tk
# ============================================================
# Assigns taxonomy to MAGs using the Genome Taxonomy Database
# (GTDB) and phylogenetic placement. Classifies genomes as
# bacteria or archaea and provides detailed taxonomic lineages.
#
# Input:  /work_meta/06_MAGs_selected/*.fa
# Output: /work_meta/08_gtdbtk/
#         - gtdbtk.bac120.summary.tsv (bacterial taxonomy)
#         - gtdbtk.ar53.summary.tsv (archaeal taxonomy)
#         - classify/ (detailed classification results)
#         - align/ (aligned sequences)
#         - identify/ (marker gene identification)
#
# Tools: GTDB-Tk v2.5.2 (classify_wf workflow)
#
# Usage: sbatch 8a_taxonomy_gtdbtk.sh
# Runtime: ~24-48 hours depending on number of MAGs
# ============================================================

# ======== BASE PATHS ========
GENOME_DIR="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected"
OUT_DIR="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/08_gtdbtk"
LOG_DIR="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs/gtdbtk"

mkdir -p "${OUT_DIR}"
mkdir -p "${LOG_DIR}"

# ======== LOAD ENVIRONMENT ========
module purge
module load mamba/24.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gtdbtk-2.5.2

echo "[$(date)] Using environment: gtdbtk-2.5.2"
echo "[$(date)] gtdbtk version:"
gtdbtk -v || { echo "[$(date)] ERROR: gtdbtk cannot be executed"; exit 1; }

echo "[$(date)] GTDBTK_DATA_PATH = ${GTDBTK_DATA_PATH:-NOT_DEFINED}"

# ======== CHECK MAGs EXIST ========
if [ ! -d "${GENOME_DIR}" ]; then
  echo "[$(date)] ERROR: GENOME_DIR does not exist: ${GENOME_DIR}"
  exit 1
fi

N_GENOMES=$(ls "${GENOME_DIR}"/*.fa 2>/dev/null | wc -l || true)

if [ "${N_GENOMES}" -eq 0 ]; then
  echo "[$(date)] ERROR: No .fa files found in ${GENOME_DIR}"
  exit 1
fi

echo "[$(date)] Found ${N_GENOMES} MAGs in ${GENOME_DIR}"
echo "[$(date)] Output will go to: ${OUT_DIR}"

# ======== RUN GTDB-TK ========

gtdbtk classify_wf \
  --genome_dir "${GENOME_DIR}" \
  --out_dir "${OUT_DIR}" \
  --cpus "${SLURM_CPUS_PER_TASK}" \
  --pplacer_cpus "${SLURM_CPUS_PER_TASK}" \
  --extension fa

echo "[$(date)] gtdbtk classify_wf finished."

# ======== CHECK OUTPUT ========
SUMMARY="${OUT_DIR}/gtdbtk.bac120.summary.tsv"
if [ -f "${SUMMARY}" ]; then
  echo "[$(date)] Summary found: ${SUMMARY}"
  head -5 "${SUMMARY}"
else
  echo "[$(date)] WARNING: ${SUMMARY} not found, check logs in ${OUT_DIR}"
fi

echo "[$(date)] ==== End of gtdbtk classify_wf ===="
