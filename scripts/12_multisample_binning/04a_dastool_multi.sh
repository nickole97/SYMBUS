#!/bin/bash
#SBATCH --job-name=ms_dastool
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/dastool_multi_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/dastool_multi_%A_%a.err

# ============================================================
# 04a_dastool_multi.sh
# Author: Nickole Villabona
# DAS Tool refinement of the MULTI-sample bins
# ============================================================
# Integrates the three multi-sample binner outputs. Identical parameters to
# 4d_refine_dastool.sh (diamond search engine, --write_bins) so the only
# difference versus the baseline arm is the depth information the binners
# were given.
#
# The per-binner provenance of every output bin is recorded by DAS Tool in
# the bin_set column of its summary table — that column is what 07 uses to
# test the CONCOCT retention hypothesis.
#
# Input:  ${RESULTS}/multi/<asm>/{metabat,maxbin,concoct/bins}
#         ${ASSEMBLY_DIR}/<asm>/final.contigs.fa (read-only)
# Output: ${RESULTS}/multi/<asm>/DASTool/<asm>_DASTool_DASTool_bins/*.fa
#         ${RESULTS}/multi/<asm>/DASTool/<asm>_DASTool_DASTool_summary.tsv
#
# Usage: sbatch --array=1-22%20 04a_dastool_multi.sh
# ============================================================

# Conda activation must come BEFORE `set -euo pipefail`: activation fails
# under SLURM otherwise. This is a known quirk on HPC3.
module purge
module load miniconda3/23.5.2
source ~/.bashrc
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate binning_env
module load samtools/1.15.1

set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
check_subset

: "${SLURM_ARRAY_TASK_ID:?Run with --array=1-${N_SUBSET}}"
ASM="$(sample_from_index "${SLURM_ARRAY_TASK_ID}")"
[[ -n "${ASM}" ]] || { echo "[$(date)] No sample at index ${SLURM_ARRAY_TASK_ID} — exiting."; exit 0; }

CONTIGS="${ASSEMBLY_DIR}/${ASM}/final.contigs.fa"
ARM="${RESULTS}/multi/${ASM}"
DT="${ARM}/DASTool"
PREFIX="${DT}/${ASM}_DASTool"

[[ -s "${CONTIGS}" ]] || die "missing assembly: ${CONTIGS}"
mkdir -p "${DT}" "${LOG_DIR}"

if [[ -s "${PREFIX}_DASTool_summary.tsv" ]]; then
    echo "[$(date)] DAS Tool already ran for ${ASM} (multi), exiting."
    exit 0
fi

cd "${DT}"

# ======== CONTIG2BIN TABLES, SKIPPING EMPTY BINNERS ========
# Several subset samples legitimately produce zero bins from a given binner,
# so absent output is not an error — but zero usable binners is.
INPUTS=(); LABELS=()

add_binner() {           # add_binner <label> <dir> <extension>
    local label="$1" dir="$2" ext="$3"
    local tsv="${DT}/${ASM}_${label}_scaffolds2bin.tsv"
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

echo "[$(date)] collecting multi-sample bins for ${ASM}"
add_binner metabat2 "${ARM}/metabat"      fa
add_binner maxbin2  "${ARM}/maxbin"       fasta
add_binner concoct  "${ARM}/concoct/bins" fa

if [[ "${#INPUTS[@]}" -lt 1 ]]; then
    echo "[$(date)] no binner produced usable bins for ${ASM} (multi) — nothing to refine."
    # Leave a marker so 07 can tell "no bins" apart from "job never ran".
    : > "${DT}/NO_BINS"
    exit 0
fi

INPUT_STR=$(IFS=,; echo "${INPUTS[*]}")
LABEL_STR=$(IFS=,; echo "${LABELS[*]}")
echo "[$(date)] DAS Tool on ${ASM} with: ${LABEL_STR}"

"${DASTOOL}" \
  -i "${INPUT_STR}" \
  -l "${LABEL_STR}" \
  -c "${CONTIGS}" \
  -o "${PREFIX}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --search_engine diamond \
  --write_bins

N_BINS=$(ls "${PREFIX}_DASTool_bins"/*.fa 2>/dev/null | wc -l || true)
echo "[$(date)] DONE ${ASM} (multi): ${N_BINS} refined bins"
[[ -s "${PREFIX}_DASTool_summary.tsv" ]] && \
    awk -F'\t' 'NR>1 {print $2}' "${PREFIX}_DASTool_summary.tsv" | sort | uniq -c
