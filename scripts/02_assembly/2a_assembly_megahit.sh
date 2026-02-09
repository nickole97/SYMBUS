#!/usr/bin/env bash
#SBATCH --job-name=assembly
#SBATCH --output=assembly.o%A_%a
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --partition=standard

# ============================================================
# 2a_assembly_megahit.sh
# Author: Nickole Villabona
# Metagenomic assembly with MEGAHIT
# ============================================================
# Assembles non-host reads into contigs using MEGAHIT assembler.
# Uses both paired and singleton reads when available.
#
# Input:  /work_meta/preprocessing/02_host_removal/<sample>/
#         - *_R1_nonhost.fastq.gz, *_R2_nonhost.fastq.gz
#         - *_singleton_nonhost.fastq.gz (optional)
# Output: /work_meta/03_assembly/<sample>/
#         - final.contigs.fa (minimum length 1000 bp)
#
# Tools: MEGAHIT v1.2.9 (min-contig-len=1000, 4 threads)
#
# Usage: sbatch --array=0-112 2a_assembly_megahit.sh
# Runtime: ~12-24 hours per sample (varies with read depth)
# ============================================================

set -euo pipefail

# ======== MODULES ========
module purge
module load megahit/1.2.9

# ======== PATHS ========
BASE="/dfs10/hammert-lab/nvillabo/Ch2/work_meta"
IN="${BASE}/preprocessing/02_host_removal"
OUT="${BASE}/03_assembly"
SAMPLES="${BASE}/preprocessing/samples.txt"

mkdir -p "${OUT}/logs"

# ======== ARRAY: get sample ID based on index ========
: "${SLURM_ARRAY_TASK_ID:?Run with --array=0-(N-1)}"
SID="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${SAMPLES}" || true)"
if [[ -z "${SID}" ]]; then
  echo "[ARRAY] No sample for index ${SLURM_ARRAY_TASK_ID}; exiting."
  exit 0
fi

# ======== INPUT FILES ========
R1="${IN}/${SID}/${SID}_R1_nonhost.fastq.gz"
R2="${IN}/${SID}/${SID}_R2_nonhost.fastq.gz"
RS="${IN}/${SID}/${SID}_singleton_nonhost.fastq.gz"
OD="${OUT}/${SID}"
LOG="${OUT}/logs/megahit_${SID}.log"

# ======== CHECKS ========
if [[ -s "${OD}/final.contigs.fa" ]]; then
  echo "✅ Assembly already exists: ${SID} → ${OD}/final.contigs.fa"
  exit 0
fi

if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
  echo "⚠️ Missing non-host paired FASTQs for ${SID}; skipping."
  echo "  R1=${R1}"
  echo "  R2=${R2}"
  exit 0
fi

# ======== ASSEMBLY ========
echo "[INFO] Assembling ${SID} | $(date)"
echo "[INFO] Threads=${SLURM_CPUS_PER_TASK}  Output=${OD}"

RS_ARG=()
[[ -s "${RS}" ]] && RS_ARG=(-r "${RS}")

megahit \
  -1 "${R1}" \
  -2 "${R2}" \
  "${RS_ARG[@]}" \
  -o "${OD}" \
  --min-contig-len 1000 \
  --num-cpu-threads "${SLURM_CPUS_PER_TASK}" \
  2>&1 | tee "${LOG}"

# ======== VERIFICATION ========
if [[ -s "${OD}/final.contigs.fa" ]]; then
  echo "[DONE] ${SID} completed successfully | $(date)"
else
  echo "❌ MEGAHIT did not generate final.contigs.fa for ${SID}"
  exit 1
fi
