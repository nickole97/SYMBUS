#!/bin/bash
#SBATCH --job-name=ms_maxbin
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/maxbin_multi_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/maxbin_multi_%A_%a.err

# ============================================================
# 03b_maxbin_multi.sh
# Author: Nickole Villabona
# MaxBin2 binning with multi-sample abundance files
# ============================================================
# Differs from 4b_binning_maxbin.sh in two ways:
#
#   1. multi-sample: one abundance file per sample, passed via -abund_list.
#   2. CORRECT FORMAT: each file is 2 columns (contig, mean depth), which is
#      what MaxBin2's abundance parser actually expects. The single-sample
#      run passed the raw jgi depth table to -abund, so MaxBin2 read column 2
#      (contigLen) as the abundance value — bug (a) in the bug ledger.
#
# Because of (2) this script cannot be compared directly against the
# existing single-sample MaxBin2 bins. That is what 03d_maxbin_single_fixed.sh
# is for: it rebuilds the single-sample arm with the format fixed, so the
# only difference left between the two arms is one sample vs many.
#
# Input:  ${SCRATCH}/02_depth/<asm>/abund_list.txt (+ abund_<sample>.txt)
#         ${ASSEMBLY_DIR}/<asm>/final.contigs.fa (read-only)
# Output: ${RESULTS}/multi/<asm>/maxbin/<asm>_maxbin.*.fasta
#
# Note: 8 h wall time because MaxBin2 scales with the number of abundance
# files and is the slowest of the three binners here.
#
# Usage: sbatch --array=1-22%20 03b_maxbin_multi.sh
# ============================================================

set -euo pipefail

# Under sbatch the script runs from SLURM's spool directory, so BASH_SOURCE does
# NOT point at this directory (first run, 2026-11-25: all 22 tasks died with
# "config.sh: No such file or directory"). Use the absolute path; override with
# SCRIPT_DIR=... if the repo ever moves.
SCRIPT_DIR="${SCRIPT_DIR:-/dfs10/hammert-lab/nvillabo/Ch2/SYMBUS/scripts/12_multisample_binning}"
source "${SCRIPT_DIR}/config.sh"
check_subset

: "${SLURM_ARRAY_TASK_ID:?Run with --array=1-${N_SUBSET}}"
ASM="$(sample_from_index "${SLURM_ARRAY_TASK_ID}")"
[[ -n "${ASM}" ]] || { echo "[$(date)] No sample at index ${SLURM_ARRAY_TASK_ID} — exiting."; exit 0; }

CONTIGS="${ASSEMBLY_DIR}/${ASM}/final.contigs.fa"
ABUND_LIST="${SCRATCH}/02_depth/${ASM}/abund_list.txt"
OUT="${RESULTS}/multi/${ASM}/maxbin"

[[ -s "${CONTIGS}" ]]    || die "missing assembly: ${CONTIGS}"
[[ -s "${ABUND_LIST}" ]] || die "missing abundance list: ${ABUND_LIST} — run 02_depth_matrix.sh first"

# Guard the format: every listed file must exist and have exactly 2 columns.
while read -r AB; do
    [[ -s "${AB}" ]] || die "abundance file missing or empty: ${AB}"
    NF_AB=$(head -1 "${AB}" | awk -F'\t' '{print NF}')
    [[ "${NF_AB}" -eq 2 ]] || die "abundance file ${AB} has ${NF_AB} columns, MaxBin2 needs 2"
done < "${ABUND_LIST}"
N_AB=$(wc -l < "${ABUND_LIST}")

mkdir -p "${OUT}" "${LOG_DIR}"

if ls "${OUT}/${ASM}_maxbin."*.fasta >/dev/null 2>&1; then
    echo "[$(date)] MaxBin2 bins already exist for ${ASM}, exiting."
    exit 0
fi

# MaxBin2 writes a pile of intermediates into the working directory.
cd "${OUT}"

echo "[$(date)] MaxBin2 on ${ASM} with ${N_AB} abundance files"
"${MAXBIN2}" \
  -contig "${CONTIGS}" \
  -abund_list "${ABUND_LIST}" \
  -out "${ASM}_maxbin" \
  -thread "${SLURM_CPUS_PER_TASK}" \
  -min_contig_length "${MIN_CONTIG_LEN}"

N_BINS=$(ls "${OUT}/${ASM}_maxbin."*.fasta 2>/dev/null | wc -l)
echo "[$(date)] DONE ${ASM}: ${N_BINS} MaxBin2 bins"
