#!/bin/bash

# ============================================================
# 6a_select_mags.sh
# Author: Nickole Villabona
# Select high-quality MAGs based on CheckM results
# ============================================================
# Copies high-quality MAGs from individual sample binning folders
# to a centralized directory for dereplication. Reads a TSV file
# with sample IDs and bin IDs to select which MAGs to keep.
#
# Input:  /work_meta/MagsToUse.tsv (columns: SAMPLE, BINID)
#         /work_meta/05_binning/<sample>/DASTool/<sample>_DASTool_DASTool_bins/*.fa
# Output: /work_meta/06_MAGs_selected/<sample>__<binid>.fa
#
# Note: The TSV file should list MAGs that meet quality thresholds:
#       - High quality: Completeness ≥90%, Contamination <5%
#       - Medium quality: Completeness ≥70%, Contamination <10%
#
# Usage: bash 6a_select_mags.sh
# ============================================================

set -euo pipefail

WORK_BASE=/dfs10/hammert-lab/nvillabo/Ch2/work_meta
BIN_BASE=${WORK_BASE}/05_binning
SELECTED_TSV=${WORK_BASE}/MagsToUse.tsv
DEST_DIR=${WORK_BASE}/06_MAGs_selected

mkdir -p "${DEST_DIR}"

echo "Reading MAG list from: ${SELECTED_TSV}"
echo "Copying selected MAGs to: ${DEST_DIR}"
echo

# Skip header line
tail -n +2 "${SELECTED_TSV}" | while read -r SAMPLE BINID; do
    if [ -z "${SAMPLE}" ] || [ -z "${BINID}" ]; then
        continue
    fi

    echo "Processing: ${SAMPLE}    ${BINID}"

    # Search candidate paths in priority order
    CANDIDATES=(
        "${BIN_BASE}/${SAMPLE}/CheckM/bins/${BINID}"
        "${BIN_BASE}/${SAMPLE}/CheckM/bins/${BINID}.fa"
        "${BIN_BASE}/${SAMPLE}/CheckM/bins/${BINID}.fasta"
        "${BIN_BASE}/${SAMPLE}/metabat/${BINID}.fa"
        "${BIN_BASE}/${SAMPLE}/maxbin/${BINID}.fasta"
        "${BIN_BASE}/${SAMPLE}/concoct/bins/${BINID}.fa"
        "${BIN_BASE}/${SAMPLE}/DASTool/${SAMPLE}_DASTool_DASTool_bins/${BINID}.fa"
    )

    FOUND=""
    for C in "${CANDIDATES[@]}"; do
        if [ -f "${C}" ]; then
            FOUND="${C}"
            break
        fi
    done

    if [ -z "${FOUND}" ]; then
        echo "⚠️  NOT FOUND: ${SAMPLE} / ${BINID}"
        continue
    fi

    # Output name: SampleID__BinId.fa
    OUT="${DEST_DIR}/${SAMPLE}__${BINID}.fa"

    cp "${FOUND}" "${OUT}"

    echo "   -> Copied: ${FOUND}  ->  ${OUT}"
done

echo
echo "=== Selected MAGs saved in: ${DEST_DIR}"
