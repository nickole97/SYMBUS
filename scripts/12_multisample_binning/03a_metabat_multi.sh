#!/bin/bash
#SBATCH --job-name=ms_metabat
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/metabat_multi_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/metabat_multi_%A_%a.err

# ============================================================
# 03a_metabat_multi.sh
# Author: Nickole Villabona
# MetaBAT2 binning with the multi-sample depth matrix
# ============================================================
# Identical to 4a_binning_metabat.sh except for the depth file: MetaBAT2
# consumes the multi-column jgi table natively, so this is the one binner
# that needs no reformatting to go multi-sample.
#
# Parameters are copied from 4a (-m 1500) and the binary is the same
# binning_env MetaBAT 2.18. Changing either would confound the comparison.
#
# Input:  ${DEPTH_DIR}/<asm>/depth.txt
#         ${ASSEMBLY_DIR}/<asm>/final.contigs.fa (read-only)
# Output: ${RESULTS}/multi/<asm>/metabat/<asm>_metabat.*.fa
#
# Usage: sbatch --array=1-22%20 03a_metabat_multi.sh
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
DEPTH="${DEPTH_DIR}/${ASM}/depth.txt"
OUT="${RESULTS}/multi/${ASM}/metabat"

[[ -s "${CONTIGS}" ]] || die "missing assembly: ${CONTIGS}"
[[ -s "${DEPTH}" ]]   || die "missing multi-sample depth: ${DEPTH} — run 02_depth_matrix.sh first"

mkdir -p "${OUT}" "${LOG_DIR}"

if ls "${OUT}/${ASM}_metabat."*.fa >/dev/null 2>&1 || [[ -e "${OUT}/NO_BINS" ]]; then
    echo "[$(date)] MetaBAT2 already ran for ${ASM} (bins or NO_BINS marker present), exiting."
    exit 0
fi

N_SAMP=$(( ( $(head -1 "${DEPTH}" | awk -F'\t' '{print NF}') - 3 ) / 2 ))
echo "[$(date)] MetaBAT2 on ${ASM} with depth from ${N_SAMP} samples"

# MetaBAT2 exits non-zero when it forms 0 bins ("0 bins (0 bases in total)
# formed."). For the failure block that is a legitimate result (5 of 22 on
# 2026-11-25, the same samples with 0 MetaBAT2 bins in the original run), not
# an error: record NO_BINS and exit 0 so DAS Tool (afterok) still runs.
RUN_LOG="${OUT}/${ASM}_metabat_run.log"
set +e
"${METABAT2}" \
  -i "${CONTIGS}" \
  -a "${DEPTH}" \
  -o "${OUT}/${ASM}_metabat" \
  -m "${MIN_CONTIG_LEN}" \
  -t "${SLURM_CPUS_PER_TASK}" 2>&1 | tee "${RUN_LOG}"
MB_RC=${PIPESTATUS[0]}
set -e
if [[ "${MB_RC}" -ne 0 ]]; then
    if grep -q "0 bins (0 bases in total) formed" "${RUN_LOG}"; then
        echo "[$(date)] DONE ${ASM}: 0 MetaBAT2 bins (NO_BINS)"
        : > "${OUT}/NO_BINS"
        exit 0
    fi
    die "MetaBAT2 exited ${MB_RC} on ${ASM} for a reason other than '0 bins formed' — see ${RUN_LOG}"
fi

# Count with find, not `ls | wc -l`: under `set -o pipefail` an `ls` with no
# match fails the whole pipeline and errexit kills the script right here —
# which is why the 5 zero-bin samples "failed" in the original run too.
N_BINS=$(find "${OUT}" -maxdepth 1 -name "${ASM}_metabat.*.fa" | wc -l)
if [[ "${N_BINS}" -eq 0 ]]; then
    echo "[$(date)] DONE ${ASM}: 0 MetaBAT2 bins (NO_BINS)"
    : > "${OUT}/NO_BINS"
    exit 0
fi
echo "[$(date)] DONE ${ASM}: ${N_BINS} MetaBAT2 bins"
