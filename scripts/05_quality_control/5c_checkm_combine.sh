#!/bin/bash

# ============================================================
# 5c_checkm_combine.sh
# Author: Nickole Villabona
# Combine CheckM quality reports from all samples
# ============================================================
# Combines individual CheckM quality tables into a single
# master table with sample identifiers. Adds a SAMPLE column
# to track which sample each MAG came from.
#
# Input:  /work_meta/05_binning/<sample>/CheckM/<sample>_checkm_qa.tsv
# Output: /work_meta/05_binning/checkm_qa_ALL.tsv
#
# Usage: bash 5c_checkm_combine.sh
# Note: Runs locally, not as SLURM job
# ============================================================

set -euo pipefail

WORK_BASE=/dfs10/hammert-lab/nvillabo/Ch2/work_meta
BIN_DIR=${WORK_BASE}/05_binning

OUTFILE=${BIN_DIR}/checkm_qa_ALL.tsv

echo "Using BIN_DIR=${BIN_DIR}"
echo "Output file: ${OUTFILE}"

HEADER_WRITTEN=0

# Clean previous output if exists
rm -f "${OUTFILE}"

# Process each sample folder
for d in "${BIN_DIR}"/*/; do
    SAMPLE=$(basename "${d}")
    QA_FILE="${d}/CheckM/${SAMPLE}_checkm_qa.tsv"

    if [ ! -f "${QA_FILE}" ]; then
        echo "[-] Not found ${QA_FILE}, skipping ${SAMPLE}"
        continue
    fi

    echo "[+] Processing ${QA_FILE}"

    # Write header only once (from first file)
    if [ ${HEADER_WRITTEN} -eq 0 ]; then
        head -n 1 "${QA_FILE}" > "${OUTFILE}"
        HEADER_WRITTEN=1
    fi

    # Add rows (without header), adding SAMPLE column at the beginning
    awk -v s="${SAMPLE}" 'NR>1 {print s "\t" $0}' "${QA_FILE}" >> "${OUTFILE}"

done

echo "Done. Combined table in: ${OUTFILE}"
