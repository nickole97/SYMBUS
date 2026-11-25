#!/bin/bash
#SBATCH --job-name=ms_checkm
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/checkm_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/checkm_%A_%a.err

# ============================================================
# 05_checkm1_both.sh
# Author: Nickole Villabona
# CheckM1 quality assessment of both new arms
# ============================================================
# Two array tasks, one per arm (1 = single_fixed, 2 = multi). Each gathers
# every DAS Tool bin from its arm into one staging directory and runs
# CheckM once over all of them, instead of 22 separate runs.
#
# CheckM1 ON PURPOSE, not CheckM2. CheckM2 is the better tool, but the
# third arm (single-BROKEN) was measured with CheckM1 from checkm_env and
# mixing versions across arms would invalidate the comparison. Those
# existing results are reused as-is, so nothing needs rerunning for it.
#
# Bins are staged as <sample>__<binid>.fa, matching the naming convention
# already used in 06_MAGs_selected/, so the sample is recoverable from the
# bin name and bins from different samples cannot collide.
#
# Input:  ${RESULTS}/<arm>/<sample>/DASTool/<sample>_DASTool_DASTool_bins/*.fa
# Output: ${RESULTS}/checkm/<arm>/                 — CheckM working tree
#         ${RESULTS}/checkm/checkm_qa_<arm>.tsv    — SAMPLE, BIN + CheckM stats
#
# NOTE: the output table gets a CORRECT header. The existing
# 05_binning/checkm_qa_ALL.tsv is written by 5c_checkm_combine.sh, which
# prepends a SAMPLE column to the data rows but not to the header, leaving
# it shifted by one (bug (c)). That file is read positionally in 07; this
# one does not need to be.
#
# Usage: sbatch --array=1-2 05_checkm1_both.sh
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
conda activate checkm_env

set -euo pipefail

# Under sbatch the script runs from SLURM's spool directory, so BASH_SOURCE does
# NOT point at this directory (first run, 2026-11-25: all 22 tasks died with
# "config.sh: No such file or directory"). Use the absolute path; override with
# SCRIPT_DIR=... if the repo ever moves.
SCRIPT_DIR="${SCRIPT_DIR:-/dfs10/hammert-lab/nvillabo/Ch2/SYMBUS/scripts/12_multisample_binning}"
source "${SCRIPT_DIR}/config.sh"
check_subset

export CHECKM_DATA_PATH="${CHECKM_DATA}"

: "${SLURM_ARRAY_TASK_ID:?Run with --array=1-2 (1=single_fixed, 2=multi)}"
case "${SLURM_ARRAY_TASK_ID}" in
    1) ARM=single_fixed ;;
    2) ARM=multi ;;
    *) die "array index must be 1 (single_fixed) or 2 (multi), got ${SLURM_ARRAY_TASK_ID}" ;;
esac

STAGE="${RESULTS}/checkm/${ARM}/bins"
CHECKM_OUT="${RESULTS}/checkm/${ARM}/out"
QA="${RESULTS}/checkm/checkm_qa_${ARM}.tsv"

mkdir -p "${STAGE}" "${CHECKM_OUT}" "${LOG_DIR}"

echo "[$(date)] arm=${ARM} | CHECKM_DATA_PATH=${CHECKM_DATA_PATH}"

# ======== STAGE BINS ========
rm -f "${STAGE}"/*.fa
N_STAGED=0
N_NOBINS=0
while read -r SID; do
    BIN_DIR="${RESULTS}/${ARM}/${SID}/DASTool/${SID}_DASTool_DASTool_bins"
    if [[ ! -d "${BIN_DIR}" ]]; then
        if [[ -f "${RESULTS}/${ARM}/${SID}/DASTool/NO_BINS" ]]; then
            N_NOBINS=$(( N_NOBINS + 1 )); continue
        fi
        die "neither bins nor a NO_BINS marker for ${SID} in arm ${ARM} — did 04 finish?"
    fi
    for B in "${BIN_DIR}"/*.fa; do
        [[ -e "${B}" ]] || continue
        cp "${B}" "${STAGE}/${SID}__$(basename "${B}")"
        N_STAGED=$(( N_STAGED + 1 ))
    done
done < <(subset_samples)

echo "[$(date)] staged ${N_STAGED} bins from ${N_SUBSET} samples (${N_NOBINS} samples had no bins)"
if [[ "${N_STAGED}" -eq 0 ]]; then
    echo "[$(date)] nothing to assess for arm ${ARM} — writing an empty table."
    printf 'SAMPLE\tBIN\n' > "${QA}"
    exit 0
fi

# ======== CHECKM ========
if [[ ! -s "${CHECKM_OUT}/storage/bin_stats_ext.tsv" ]]; then
    echo "[$(date)] running checkm lineage_wf over ${N_STAGED} bins"
    "${CHECKM}" lineage_wf \
      -x fa \
      -t "${SLURM_CPUS_PER_TASK}" \
      --pplacer_threads "${SLURM_CPUS_PER_TASK}" \
      "${STAGE}" \
      "${CHECKM_OUT}"
else
    echo "[$(date)] CheckM output already exists, skipping lineage_wf"
fi

# ======== QA TABLE WITH A CORRECT HEADER ========
RAW="${RESULTS}/checkm/${ARM}/qa_raw.tsv"
"${CHECKM}" qa \
  -o 2 --tab_table \
  -t "${SLURM_CPUS_PER_TASK}" \
  -f "${RAW}" \
  "${CHECKM_OUT}/lineage.ms" \
  "${CHECKM_OUT}"

# Split the staged "<sample>__<bin>" name back into two real columns, and
# keep every CheckM column with its own name intact.
awk -F'\t' -v OFS='\t' '
  NR==1 { printf "SAMPLE\tBIN"; for (i=2; i<=NF; i++) printf "%s%s", OFS, $i; printf "\n"; next }
  {
    split($1, p, "__")
    sample = p[1]
    bin = substr($1, length(sample) + 3)
    printf "%s\t%s", sample, bin
    for (i=2; i<=NF; i++) printf "%s%s", OFS, $i
    printf "\n"
  }
' "${RAW}" > "${QA}.tmp"
mv "${QA}.tmp" "${QA}"

echo "[$(date)] DONE ${ARM}: $(( $(wc -l < "${QA}") - 1 )) bins assessed -> ${QA}"
echo "[$(date)] MIMAG tiers (HQ >90/<5, MQ >=50/<10):"
awk -F'\t' '
  NR==1 { for (i=1;i<=NF;i++) { if ($i=="Completeness") c=i; if ($i=="Contamination") x=i } next }
  { if ($c>90 && $x<5) hq++; else if ($c>=50 && $x<10) mq++; else lq++ }
  END { printf "           HQ=%d MQ=%d LQ=%d\n", hq+0, mq+0, lq+0 }
' "${QA}"
