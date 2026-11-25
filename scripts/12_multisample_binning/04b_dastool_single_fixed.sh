#!/bin/bash
#SBATCH --job-name=sf_dastool
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/dastool_sf_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/dastool_sf_%A_%a.err

# ============================================================
# 04b_dastool_single_fixed.sh
# Author: Nickole Villabona
# DAS Tool refinement of the BASELINE arm (single-sample, MaxBin2 fixed)
# ============================================================
# Mixes two sources on purpose:
#   MetaBAT2 + CONCOCT : the EXISTING bins from 05_binning, read-only. They
#                        were given correct input in the original run, so
#                        rerunning them would only add noise and compute.
#   MaxBin2            : the rebuilt bins from 03d, with the abundance
#                        format fixed.
#
# Nothing here writes into 05_binning. The single-sample side stays frozen;
# only the refined output goes to ${RESULTS}/single_fixed/.
#
# Input:  ${SINGLE_BIN_DIR}/<s>/metabat/*.fa          (read-only)
#         ${SINGLE_BIN_DIR}/<s>/concoct/bins/*.fa     (read-only)
#         ${RESULTS}/single_fixed/<s>/maxbin/*.fasta  (from 03d)
#         ${ASSEMBLY_DIR}/<s>/final.contigs.fa        (read-only)
# Output: ${RESULTS}/single_fixed/<s>/DASTool/<s>_DASTool_DASTool_bins/*.fa
#         ${RESULTS}/single_fixed/<s>/DASTool/<s>_DASTool_DASTool_summary.tsv
#
# Usage: sbatch --array=1-22%20 04b_dastool_single_fixed.sh
# ============================================================

# Conda before `set -euo pipefail` — see 04a.
module purge
module load miniconda3/23.5.2
source ~/.bashrc
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate binning_env
module load samtools/1.15.1

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
ARM="${RESULTS}/single_fixed/${SID}"
DT="${ARM}/DASTool"
PREFIX="${DT}/${SID}_DASTool"

[[ -s "${CONTIGS}" ]] || die "missing assembly: ${CONTIGS}"

# 03d must have run, even if it produced no bins — otherwise this arm would
# silently be "MetaBAT2 + CONCOCT only" and not comparable to the multi arm.
[[ -d "${ARM}/maxbin" ]] || die "no single_fixed/maxbin directory for ${SID} — run 03d_maxbin_single_fixed.sh first"

mkdir -p "${DT}" "${LOG_DIR}"

if [[ -s "${PREFIX}_DASTool_summary.tsv" ]]; then
    echo "[$(date)] DAS Tool already ran for ${SID} (single-fixed), exiting."
    exit 0
fi

cd "${DT}"

INPUTS=(); LABELS=()

add_binner() {           # add_binner <label> <dir> <extension>
    local label="$1" dir="$2" ext="$3"
    local tsv="${DT}/${SID}_${label}_scaffolds2bin.tsv"
    if ls "${dir}"/*."${ext}" >/dev/null 2>&1; then
        "${FASTA2BIN}" -e "${ext}" -i "${dir}" > "${tsv}"
        if [[ -s "${tsv}" ]]; then
            INPUTS+=("${tsv}"); LABELS+=("${label}")
            echo "[$(date)]   ${label}: $(cut -f2 "${tsv}" | sort -u | wc -l) bins"
        else
            echo "[$(date)]   ${label}: contig2bin table empty, skipping"
        fi
    else
        echo "[$(date)]   ${label}: no bins, skipping"
    fi
}

echo "[$(date)] collecting single-sample bins for ${SID}"
echo "[$(date)]   (metabat2 + concoct reused from ${SINGLE_BIN_DIR}, read-only)"
add_binner metabat2 "${SINGLE_BIN_DIR}/${SID}/metabat"      fa
add_binner maxbin2  "${ARM}/maxbin"                          fasta
add_binner concoct  "${SINGLE_BIN_DIR}/${SID}/concoct/bins"  fa

if [[ "${#INPUTS[@]}" -lt 1 ]]; then
    echo "[$(date)] no binner produced usable bins for ${SID} (single-fixed) — nothing to refine."
    : > "${DT}/NO_BINS"
    exit 0
fi

INPUT_STR=$(IFS=,; echo "${INPUTS[*]}")
LABEL_STR=$(IFS=,; echo "${LABELS[*]}")
echo "[$(date)] DAS Tool on ${SID} with: ${LABEL_STR}"

"${DASTOOL}" \
  -i "${INPUT_STR}" \
  -l "${LABEL_STR}" \
  -c "${CONTIGS}" \
  -o "${PREFIX}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --search_engine diamond \
  --write_bins

N_BINS=$(ls "${PREFIX}_DASTool_bins"/*.fa 2>/dev/null | wc -l || true)
echo "[$(date)] DONE ${SID} (single-fixed): ${N_BINS} refined bins"
[[ -s "${PREFIX}_DASTool_summary.tsv" ]] && \
    awk -F'\t' 'NR>1 {print $2}' "${PREFIX}_DASTool_summary.tsv" | sort | uniq -c
