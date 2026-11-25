#!/bin/bash
#SBATCH --job-name=ms_concoct
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/concoct_multi_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/concoct_multi_%A_%a.err

# ============================================================
# 03c_concoct_multi.sh
# Author: Nickole Villabona
# CONCOCT binning with the multi-sample coverage table
# ============================================================
# This is the headline test of the whole experiment. CONCOCT is designed to
# use coverage covariance ACROSS samples, but the single-sample pipeline gave
# it exactly one coverage column. It is still the largest contributor to the
# final MAG set (86 of 183 DAS Tool bins, 59 of them high quality), so the
# hypothesis is not that it starts working — it is that its RETENTION RATE
# improves: it currently produces 16-55 bins per sample of which DAS Tool
# keeps 0-2.
#
# The cut-up FASTA and the coverage table are both built in 02 (the only
# step allowed to read the BAMs), so this script just clusters.
#
# Input:  ${DEPTH_DIR}/<asm>/contigs_10K.fa
#         ${DEPTH_DIR}/<asm>/coverage_table.tsv
#         ${ASSEMBLY_DIR}/<asm>/final.contigs.fa (read-only, for bin extraction)
# Output: ${RESULTS}/multi/<asm>/concoct/bins/*.fa
#
# Note: 12 h wall time — CONCOCT's VBGMM scales with chunks x samples, and
# here it gets 22 coverage columns instead of 1.
#
# Usage: sbatch --array=1-22%20 03c_concoct_multi.sh
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

module load samtools/1.15.1

: "${SLURM_ARRAY_TASK_ID:?Run with --array=1-${N_SUBSET}}"
ASM="$(sample_from_index "${SLURM_ARRAY_TASK_ID}")"
[[ -n "${ASM}" ]] || { echo "[$(date)] No sample at index ${SLURM_ARRAY_TASK_ID} — exiting."; exit 0; }

CONTIGS="${ASSEMBLY_DIR}/${ASM}/final.contigs.fa"
DEPTH_DIR="${DEPTH_DIR}/${ASM}"
CUT_FA="${DEPTH_DIR}/contigs_10K.fa"
COV="${DEPTH_DIR}/coverage_table.tsv"
OUT="${RESULTS}/multi/${ASM}/concoct"
BINS="${OUT}/bins"

[[ -s "${CONTIGS}" ]] || die "missing assembly: ${CONTIGS}"
[[ -s "${CUT_FA}" ]]  || die "missing cut-up FASTA: ${CUT_FA} — run 02_depth_matrix.sh first"
[[ -s "${COV}" ]]     || die "missing coverage table: ${COV} — run 02_depth_matrix.sh first"

mkdir -p "${OUT}" "${LOG_DIR}"

if ls "${BINS}"/*.fa >/dev/null 2>&1; then
    echo "[$(date)] CONCOCT bins already exist for ${ASM}, exiting."
    exit 0
fi

cd "${OUT}"

N_SAMP=$(( $(head -1 "${COV}" | awk -F'\t' '{print NF}') - 1 ))
echo "[$(date)] CONCOCT on ${ASM}: $(( $(wc -l < "${COV}") - 1 )) chunks x ${N_SAMP} samples"

# ======== 1) CLUSTER ========
PREFIX=concoct_out
CLUSTER_RAW="${PREFIX}_clustering_gt1000.csv"
if [[ ! -s "${CLUSTER_RAW}" ]]; then
    "${CONCOCT}" \
      --composition_file "${CUT_FA}" \
      --coverage_file "${COV}" \
      --threads "${SLURM_CPUS_PER_TASK}" \
      -b "${PREFIX}"
else
    echo "[$(date)] raw clustering exists, skipping CONCOCT"
fi
[[ -s "${CLUSTER_RAW}" ]] || die "CONCOCT produced no clustering (${CLUSTER_RAW})"

# ======== 2) MERGE CHUNKS BACK TO WHOLE CONTIGS ========
MERGED="${PREFIX}_clustering_merged.csv"
if [[ ! -s "${MERGED}" ]]; then
    "${CONCOCT_MERGE}" "${CLUSTER_RAW}" > "${MERGED}.tmp"
    mv "${MERGED}.tmp" "${MERGED}"
else
    echo "[$(date)] merged clustering exists, skipping"
fi
[[ -s "${MERGED}" ]] || die "merge_cutup_clustering produced nothing"

# ======== 3) EXTRACT BINS ========
mkdir -p "${BINS}"
"${CONCOCT_EXTRACT}" "${CONTIGS}" "${MERGED}" --output_path "${BINS}"

N_BINS=$(find "${OUT}/bins" -maxdepth 1 -name "*.fa" 2>/dev/null | wc -l)   # find, not ls|wc: pipefail-safe
echo "[$(date)] DONE ${ASM}: ${N_BINS} CONCOCT bins (DAS Tool will keep a subset)"
