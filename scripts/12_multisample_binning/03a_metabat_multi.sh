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
# Input:  ${SCRATCH}/02_depth/<asm>/depth.txt
#         ${ASSEMBLY_DIR}/<asm>/final.contigs.fa (read-only)
# Output: ${RESULTS}/multi/<asm>/metabat/<asm>_metabat.*.fa
#
# Usage: sbatch --array=1-22%20 03a_metabat_multi.sh
# ============================================================

set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
check_subset

: "${SLURM_ARRAY_TASK_ID:?Run with --array=1-${N_SUBSET}}"
ASM="$(sample_from_index "${SLURM_ARRAY_TASK_ID}")"
[[ -n "${ASM}" ]] || { echo "[$(date)] No sample at index ${SLURM_ARRAY_TASK_ID} — exiting."; exit 0; }

CONTIGS="${ASSEMBLY_DIR}/${ASM}/final.contigs.fa"
DEPTH="${SCRATCH}/02_depth/${ASM}/depth.txt"
OUT="${RESULTS}/multi/${ASM}/metabat"

[[ -s "${CONTIGS}" ]] || die "missing assembly: ${CONTIGS}"
[[ -s "${DEPTH}" ]]   || die "missing multi-sample depth: ${DEPTH} — run 02_depth_matrix.sh first"

mkdir -p "${OUT}" "${LOG_DIR}"

if ls "${OUT}/${ASM}_metabat."*.fa >/dev/null 2>&1; then
    echo "[$(date)] MetaBAT2 bins already exist for ${ASM}, exiting."
    exit 0
fi

N_SAMP=$(( ( $(head -1 "${DEPTH}" | awk -F'\t' '{print NF}') - 3 ) / 2 ))
echo "[$(date)] MetaBAT2 on ${ASM} with depth from ${N_SAMP} samples"

"${METABAT2}" \
  -i "${CONTIGS}" \
  -a "${DEPTH}" \
  -o "${OUT}/${ASM}_metabat" \
  -m "${MIN_CONTIG_LEN}" \
  -t "${SLURM_CPUS_PER_TASK}"

N_BINS=$(ls "${OUT}/${ASM}_metabat."*.fa 2>/dev/null | wc -l)
echo "[$(date)] DONE ${ASM}: ${N_BINS} MetaBAT2 bins"
