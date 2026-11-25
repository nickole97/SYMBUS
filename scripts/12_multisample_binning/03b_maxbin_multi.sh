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
# Input:  ${DEPTH_DIR}/<asm>/abund_list.txt (+ abund_<sample>.txt)
#         ${ASSEMBLY_DIR}/<asm>/final.contigs.fa (read-only)
# Output: ${RESULTS}/multi/<asm>/maxbin/<asm>_maxbin.*.fasta
#
# Note: 8 h wall time because MaxBin2 scales with the number of abundance
# files and is the slowest of the three binners here.
#
# Usage: sbatch --array=1-22%20 03b_maxbin_multi.sh
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
ASM="$(sample_from_index "${SLURM_ARRAY_TASK_ID}")"
[[ -n "${ASM}" ]] || { echo "[$(date)] No sample at index ${SLURM_ARRAY_TASK_ID} — exiting."; exit 0; }

CONTIGS="${ASSEMBLY_DIR}/${ASM}/final.contigs.fa"
ABUND_LIST="${DEPTH_DIR}/${ASM}/abund_list.txt"
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

if ls "${OUT}/${ASM}_maxbin."*.fasta >/dev/null 2>&1 || [[ -e "${OUT}/NO_BINS" ]]; then
    echo "[$(date)] MaxBin2 already ran for ${ASM} (bins or NO_BINS marker present), exiting."
    exit 0
fi

# MaxBin2 writes a pile of intermediates into the working directory.
cd "${OUT}"

echo "[$(date)] MaxBin2 on ${ASM} with ${N_AB} abundance files"
# MaxBin2 exits 255 when the assembly has too few marker genes to seed bins
# ("the medium of marker gene number <= 1"). In the single-FIXED arm this hit
# 12 of 22 samples (the failure block). It is a legitimate zero-bin result:
# record it and exit 0 so DAS Tool (afterok) still runs; this is also one of
# the things the experiment measures (does multi-sample depth change it?).
set +e
"${MAXBIN2}" \
  -contig "${CONTIGS}" \
  -abund_list "${ABUND_LIST}" \
  -out "${ASM}_maxbin" \
  -thread "${SLURM_CPUS_PER_TASK}" \
  -min_contig_length "${MIN_CONTIG_LEN}"
MAXBIN_RC=$?
set -e
if [[ "${MAXBIN_RC}" -ne 0 ]]; then
    if grep -q "cannot be binned" "${OUT}/${ASM}_maxbin.log" 2>/dev/null; then
        echo "[$(date)] DONE ${ASM}: 0 bins — MaxBin2 found too few marker genes to seed bins (NO_BINS)"
        : > "${OUT}/NO_BINS"
        exit 0
    fi
    die "MaxBin2 exited ${MAXBIN_RC} on ${ASM} for a reason other than 'cannot be binned' — see ${OUT}/${ASM}_maxbin.log"
fi

N_BINS=$(ls "${OUT}/${ASM}_maxbin."*.fasta 2>/dev/null | wc -l)
echo "[$(date)] DONE ${ASM}: ${N_BINS} MaxBin2 bins"
