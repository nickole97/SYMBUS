#!/bin/bash
#SBATCH --job-name=ms_gtdbtk
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=08:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/gtdbtk_%j.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/gtdbtk_%j.err

# ============================================================
# 06_gtdbtk_both.sh
# Author: Nickole Villabona
# GTDB-Tk classification of both new arms in a single run
# ============================================================
# This produces THE metric the experiment exists to answer: does
# multi-sample binning recover taxa that single-sample binning missed, or
# just more copies of the same ones?
#
# Both arms go into ONE classify_wf run, with bins prefixed sf__ / ms__,
# because ~15 of the ~16 min runtime of a GTDB-Tk run on this reference set
# is a fixed cost (sketching 143k reference genomes). Two separate runs
# would pay it twice for no benefit.
#
# The existing single-BROKEN bins for these samples are already classified
# in work_meta/08_gtdbtk/gtdbtk.bac120.summary.tsv (same GTDB r226, same
# GTDB-Tk 2.5.2), so that third arm needs no run here.
#
# Input:  ${RESULTS}/{single_fixed,multi}/<s>/DASTool/<s>_DASTool_DASTool_bins/*.fa
# Output: ${RESULTS}/gtdbtk/gtdbtk.bac120.summary.tsv  (+ ar53 if any archaea)
#         ${RESULTS}/gtdbtk/genomes/{sf,ms}__<sample>__<bin>.fa
#
# MEMORY: 256 GB is deliberate and should not be trimmed. A validation run
# (job 59407025) peaked at only ~110 GiB, but both its genomes were resolved
# by the skani ANI pre-screen, which skips identify/align/pplacer entirely.
# Novel bins — exactly what this experiment hopes to recover — will not hit
# that short-circuit and will take the pplacer path, which is the expensive
# one. Same reasoning as the note in 08_annotation/README.md for 8a.
#
# Usage: sbatch 06_gtdbtk_both.sh
# ============================================================

# Conda before `set -euo pipefail`.
module purge
module load mamba/24.3.0
source ~/.bashrc
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${GTDBTK_ENV:-gtdbtk-2.5.2}"

set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
check_subset

GENOMES="${RESULTS}/gtdbtk/genomes"
OUT="${RESULTS}/gtdbtk/out"
mkdir -p "${GENOMES}" "${LOG_DIR}"

echo "[$(date)] GTDB-Tk $(gtdbtk --version 2>&1 | head -1)"
echo "[$(date)] GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH:-NOT_DEFINED}"
[[ -n "${GTDBTK_DATA_PATH:-}" ]] || die "GTDBTK_DATA_PATH is not set by the env"

# ======== STAGE BINS FROM BOTH ARMS ========
rm -f "${GENOMES}"/*.fa
declare -A STAGED=( [single_fixed]=0 [multi]=0 )

for ARM in single_fixed multi; do
    case "${ARM}" in
        single_fixed) TAG=sf ;;
        multi)        TAG=ms ;;
    esac
    while read -r SID; do
        BIN_DIR="${RESULTS}/${ARM}/${SID}/DASTool/${SID}_DASTool_DASTool_bins"
        [[ -d "${BIN_DIR}" ]] || continue
        for B in "${BIN_DIR}"/*.fa; do
            [[ -e "${B}" ]] || continue
            cp "${B}" "${GENOMES}/${TAG}__${SID}__$(basename "${B}")"
            STAGED[${ARM}]=$(( STAGED[${ARM}] + 1 ))
        done
    done < <(subset_samples)
done

N_TOTAL=$(ls "${GENOMES}"/*.fa 2>/dev/null | wc -l)
echo "[$(date)] staged ${N_TOTAL} bins: single_fixed=${STAGED[single_fixed]} multi=${STAGED[multi]}"
[[ "${N_TOTAL}" -gt 0 ]] || die "no bins staged — did 04a/04b finish?"

# ======== CLASSIFY ========
# Same invocation as 8a_taxonomy_gtdbtk.sh, so classifications are directly
# comparable to the existing single-BROKEN results.
gtdbtk classify_wf \
  --genome_dir "${GENOMES}" \
  --out_dir "${OUT}" \
  --cpus "${SLURM_CPUS_PER_TASK}" \
  --pplacer_cpus "${SLURM_CPUS_PER_TASK}" \
  --extension fa

# ======== COLLECT ========
for LEVEL in bac120 ar53; do
    SRC="${OUT}/gtdbtk.${LEVEL}.summary.tsv"
    if [[ -s "${SRC}" ]]; then
        cp "${SRC}" "${RESULTS}/gtdbtk/gtdbtk.${LEVEL}.summary.tsv"
        echo "[$(date)] ${LEVEL}: $(( $(wc -l < "${SRC}") - 1 )) genomes classified"
    fi
done

SUMMARY="${RESULTS}/gtdbtk/gtdbtk.bac120.summary.tsv"
[[ -s "${SUMMARY}" ]] || die "no bac120 summary produced — check the logs in ${OUT}"

echo
echo "[$(date)] distinct species per arm (quick look; 07 does this properly):"
awk -F'\t' 'NR>1 {
    split($1, a, "__"); arm = a[1]
    n = split($2, t, ";"); sp = t[n]
    if (sp != "s__") key = arm "\t" sp; else key = arm "\t(novel/unassigned)"
    seen[key] = 1
} END { for (k in seen) print k }' "${SUMMARY}" | sort | awk -F'\t' '{c[$1]++} END {for (a in c) printf "  %s: %d distinct species-level labels\n", a, c[a]}'

echo "[$(date)] DONE"
