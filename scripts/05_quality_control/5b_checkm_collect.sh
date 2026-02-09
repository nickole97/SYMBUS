#!/bin/bash

# ============================================================
# 5b_checkm_collect.sh
# Author: Nickole Villabona
# Collect CheckM quality reports for each sample
# ============================================================
# Generates summary tables from CheckM output for each sample.
# Runs CheckM qa to create tab-delimited quality reports that
# are easier to parse and combine across samples.
#
# Input:  /work_meta/05_binning/<sample>/CheckM/lineage.ms
#         /work_meta/05_binning/<sample>/CheckM/bins/
# Output: /work_meta/05_binning/<sample>/CheckM/<sample>_checkm_qa.tsv
#
# Tools: CheckM qa (quality assessment output format 2)
#
# Usage: bash 5b_checkm_collect.sh
# Note: Runs locally, not as SLURM job
# ============================================================

WORK_BASE=/dfs10/hammert-lab/nvillabo/Ch2/work_meta
BIN_DIR=${WORK_BASE}/05_binning

# Absolute path to checkm
CHECKM=~/.conda/envs/checkm_env/bin/checkm

for SAMPLE in ${BIN_DIR}/*; do
    SAMPLE_NAME=$(basename "${SAMPLE}")
    CHECKM_DIR=${SAMPLE}/CheckM

    # Only process if CheckM folder exists
    if [ ! -d "${CHECKM_DIR}" ]; then
        continue
    fi

    echo "Processing ${SAMPLE_NAME}..."

    cd "${CHECKM_DIR}"

    # Skip if summary already exists
    if [ -f "${SAMPLE_NAME}_checkm_qa.tsv" ]; then
        echo "  ${SAMPLE_NAME}_checkm_qa.tsv already exists, skipping."
        cd - >/dev/null
        continue
    fi

    # Generate quality assessment table
    # 'bins' is the folder with MAGs that were passed to CheckM
    "$CHECKM" qa \
      lineage.ms bins \
      -o 2 \
      --tab_table \
      -t 8 \
      -f "${SAMPLE_NAME}_checkm_qa.tsv"

    cd - >/dev/null
done

echo "[DONE] CheckM QA tables generated for all samples"
