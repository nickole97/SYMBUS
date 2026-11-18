#!/bin/bash
# ============================================================
# 00_make_backmap_pairs.sh
# Author: Nickole Villabona
# Build the backmap pair list that drives the whole experiment
# ============================================================
# Calls the selector named in config.sh and writes its output to
# backmap_pairs.tsv, then validates that every referenced assembly and
# read file actually exists before any compute is spent.
#
# Input:  subset_22.tsv (via config.sh)
#         selectors/${SELECTOR}.sh
# Output: ${RESULTS}/backmap_pairs.tsv  — <assembly_sample>\t<read_sample>
#
# This file is the pluggable seam: swapping all_vs_all for simka changes
# only which pairs appear here. 01_crossmap.sh indexes it by SLURM array
# task, so it must not be regenerated once mapping has started.
#
# Usage: bash 00_make_backmap_pairs.sh
# Note: runs locally on the login node — no computation, just file checks.
# ============================================================

set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
check_subset

SELECTOR_SH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/selectors/${SELECTOR}.sh"
[[ -f "${SELECTOR_SH}" ]] || die "selector not found: ${SELECTOR_SH}"

mkdir -p "${RESULTS}" "${LOG_DIR}"

echo "[$(date)] subset:   ${N_SUBSET} samples (${SUBSET_TSV})"
echo "[$(date)] selector: ${SELECTOR}"

# ======== REFUSE TO CLOBBER A LIST ALREADY IN USE ========
# Array indices point into this file. Regenerating it after 01_crossmap.sh
# has run would silently remap task numbers to different pairs.
if [[ -s "${PAIRS_FILE}" ]]; then
    N_OLD=$(wc -l < "${PAIRS_FILE}")
    echo "[$(date)] WARNING: ${PAIRS_FILE} already exists (${N_OLD} pairs)."
    if [[ -d "${SCRATCH}/01_crossmap" ]] && \
       [[ -n "$(ls -A "${SCRATCH}/01_crossmap" 2>/dev/null)" ]]; then
        die "cross-mapping has already started against the existing pair list.
       Refusing to overwrite it — array task indices would shift.
       To start over, move ${SCRATCH}/01_crossmap and ${PAIRS_FILE} aside first."
    fi
    cp -v "${PAIRS_FILE}" "${PAIRS_FILE}.$(date +%Y%m%d%H%M%S).bak"
fi

# ======== GENERATE ========
TMP_PAIRS=$(mktemp)
trap 'rm -f "${TMP_PAIRS}"' EXIT
bash "${SELECTOR_SH}" > "${TMP_PAIRS}"

N_PAIRS=$(wc -l < "${TMP_PAIRS}")
[[ "${N_PAIRS}" -gt 0 ]] || die "selector produced no pairs"

# ======== VALIDATE ========
# Fail now, loudly, rather than 300 array tasks in. Check each distinct
# assembly and each distinct read set once (44 checks, not 484).
MISSING=0
while read -r ASM; do
    CONTIGS="${ASSEMBLY_DIR}/${ASM}/final.contigs.fa"
    [[ -s "${CONTIGS}" ]] || { echo "[MISSING] assembly: ${CONTIGS}"; MISSING=$((MISSING+1)); }
done < <(cut -f1 "${TMP_PAIRS}" | sort -u)

while read -r RD; do
    R1="${READS_DIR}/${RD}/${RD}_R1_nonhost.fastq.gz"
    R2="${READS_DIR}/${RD}/${RD}_R2_nonhost.fastq.gz"
    [[ -s "${R1}" ]] || { echo "[MISSING] reads R1: ${R1}"; MISSING=$((MISSING+1)); }
    [[ -s "${R2}" ]] || { echo "[MISSING] reads R2: ${R2}"; MISSING=$((MISSING+1)); }
done < <(cut -f2 "${TMP_PAIRS}" | sort -u)

[[ "${MISSING}" -eq 0 ]] || die "${MISSING} missing input file(s) — fix before submitting 01_crossmap.sh"

# ======== WRITE ========
mv "${TMP_PAIRS}" "${PAIRS_FILE}"
trap - EXIT

N_ASM=$(cut -f1 "${PAIRS_FILE}" | sort -u | wc -l)
echo "[$(date)] wrote ${N_PAIRS} pairs over ${N_ASM} assemblies -> ${PAIRS_FILE}"
echo
echo "Next:"
echo "  sbatch --array=1-${N_PAIRS}%${THROTTLE} 01_crossmap.sh"
