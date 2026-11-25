#!/bin/bash
#SBATCH --job-name=sf_maxbin
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/maxbin_sf_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/maxbin_sf_%A_%a.err

# ============================================================
# 03d_maxbin_single_fixed.sh
# Author: Nickole Villabona
# MaxBin2 on SINGLE-sample depth, with the abundance format fixed
# ============================================================
# Builds the BASELINE arm. MaxBin2 contributes 18.6% of the existing DAS
# Tool bins (34 of 183, 16 of them high quality, and it is the only source
# of bins in 8 samples), so the single-sample arm cannot be left with the
# broken abundance format — otherwise "multi beats single" would be partly
# explained by the bug fix rather than by the method.
#
# This is the ONLY binner rerun for the single-sample arm. MetaBAT2 and
# CONCOCT received correct input in the original run, so 04b reuses their
# existing bins straight from 05_binning (read-only).
#
# The fix: the existing /04_mapping/<s>/<s>_depth.txt is a 5-column jgi
# table (contigName, contigLen, totalAvgDepth, <bam>, <bam>-var). Column 4
# is the mean depth. 4b_binning_maxbin.sh passed this whole file to -abund,
# which expects 2 columns, so MaxBin2 read contigLen as the abundance.
# Here we extract (contigName, column 4) into a proper 2-column file.
#
# Input:  ${SINGLE_MAP_DIR}/<s>/<s>_depth.txt  (read-only, NOT modified)
#         ${ASSEMBLY_DIR}/<s>/final.contigs.fa (read-only)
# Output: ${RESULTS}/single_fixed/<s>/maxbin/<s>_maxbin.*.fasta
#         ${RESULTS}/single_fixed/<s>/maxbin/abund_<s>.txt
#
# Usage: sbatch --array=1-22%20 03d_maxbin_single_fixed.sh
# ============================================================

# Conda activation must come BEFORE `set -euo pipefail` (see 04a). Tools are
# called by absolute path, but MaxBin2 / CONCOCT / CheckM need their env on
# PATH to find their own dependencies (first run 2026-11-25 died with
# "Cannot run HMMER3"). The original 4b only worked because it was submitted
# from a shell that already had the env active and sbatch inherits PATH.
module purge
module load miniconda3/23.5.2
source ~/.bashrc
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate binning_env

set -euo pipefail

# Under sbatch the script runs from SLURM's spool directory, so BASH_SOURCE does
# NOT point at this directory (first run, 2026-11-25: all 22 tasks died with
# "config.sh: No such file or directory"). Use the absolute path; override with
# SCRIPT_DIR=... if the repo ever moves.
SCRIPT_DIR="${SCRIPT_DIR:-/dfs10/hammert-lab/nvillabo/Ch2/SYMBUS/scripts/12_multisample_binning}"
source "${SCRIPT_DIR}/config.sh"
check_subset

: "${SLURM_ARRAY_TASK_ID:?Run with --array=1-${N_SUBSET}}"
SID="$(sample_from_index "${SLURM_ARRAY_TASK_ID}")"
[[ -n "${SID}" ]] || { echo "[$(date)] No sample at index ${SLURM_ARRAY_TASK_ID} — exiting."; exit 0; }

CONTIGS="${ASSEMBLY_DIR}/${SID}/final.contigs.fa"
SINGLE_DEPTH="${SINGLE_MAP_DIR}/${SID}/${SID}_depth.txt"
OUT="${RESULTS}/single_fixed/${SID}/maxbin"
ABUND="${OUT}/abund_${SID}.txt"

[[ -s "${CONTIGS}" ]]      || die "missing assembly: ${CONTIGS}"
[[ -s "${SINGLE_DEPTH}" ]] || die "missing single-sample depth: ${SINGLE_DEPTH}"

mkdir -p "${OUT}" "${LOG_DIR}"

if ls "${OUT}/${SID}_maxbin."*.fasta >/dev/null 2>&1 || [[ -e "${OUT}/NO_BINS" ]]; then
    echo "[$(date)] fixed MaxBin2 already ran for ${SID} (bins or NO_BINS marker present), exiting."
    exit 0
fi

# ======== DERIVE A CORRECT 2-COLUMN ABUNDANCE FILE ========
# Assert the layout before trusting column 4. A single-BAM jgi table must
# have exactly 5 columns; anything else means this file was not produced the
# way 3a produced it, and we must not guess.
NF_DEPTH=$(head -1 "${SINGLE_DEPTH}" | awk -F'\t' '{print NF}')
[[ "${NF_DEPTH}" -eq 5 ]] || \
    die "${SINGLE_DEPTH} has ${NF_DEPTH} columns, expected 5 (single-BAM jgi table). Refusing to guess which column is the mean depth."

awk -F'\t' 'NR>1 {printf "%s\t%s\n", $1, $4}' "${SINGLE_DEPTH}" > "${ABUND}.tmp"
mv "${ABUND}.tmp" "${ABUND}"
[[ -s "${ABUND}" ]] || die "derived abundance file is empty: ${ABUND}"

# Sanity check the fix actually did something: if column 4 equalled column 2
# the file would be identical to the broken interpretation.
IDENTICAL=$(awk -F'\t' 'NR>1 && $2==$4 {n++} END {print n+0}' "${SINGLE_DEPTH}")
N_ROWS=$(( $(wc -l < "${SINGLE_DEPTH}") - 1 ))
echo "[$(date)] abundance derived from column 4: ${N_ROWS} contigs (${IDENTICAL} where contigLen == depth, expected ~0)"

# ======== RUN MAXBIN2 ========
cd "${OUT}"

echo "[$(date)] MaxBin2 on ${SID}, single-sample depth, CORRECT 2-column abundance"
# MaxBin2 exits 255 when the assembly has too few marker genes to seed bins
# ("the medium of marker gene number <= 1"). That is a legitimate zero-bin
# result for the failure block (12 of 22 samples on 2026-11-25), not an error:
# record it and exit 0 so the DAS Tool step (afterok) can still run. Any other
# non-zero exit is a real failure.
set +e
"${MAXBIN2}" \
  -contig "${CONTIGS}" \
  -abund "${ABUND}" \
  -out "${SID}_maxbin" \
  -thread "${SLURM_CPUS_PER_TASK}" \
  -min_contig_length "${MIN_CONTIG_LEN}"
MAXBIN_RC=$?
set -e
if [[ "${MAXBIN_RC}" -ne 0 ]]; then
    if grep -q "cannot be binned" "${OUT}/${SID}_maxbin.log" 2>/dev/null; then
        echo "[$(date)] DONE ${SID}: 0 bins — MaxBin2 found too few marker genes to seed bins (NO_BINS)"
        : > "${OUT}/NO_BINS"
        exit 0
    fi
    die "MaxBin2 exited ${MAXBIN_RC} on ${SID} for a reason other than 'cannot be binned' — see ${OUT}/${SID}_maxbin.log"
fi

N_BINS=$(ls "${OUT}/${SID}_maxbin."*.fasta 2>/dev/null | wc -l)
N_OLD=$(ls "${SINGLE_BIN_DIR}/${SID}/maxbin/${SID}_maxbin."*.fasta 2>/dev/null | wc -l || true)
echo "[$(date)] DONE ${SID}: ${N_BINS} bins with fixed abundance (broken run had ${N_OLD})"
