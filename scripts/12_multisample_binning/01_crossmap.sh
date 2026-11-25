#!/bin/bash
#SBATCH --job-name=ms_crossmap
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/crossmap_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/crossmap_%A_%a.err

# ============================================================
# 01_crossmap.sh
# Author: Nickole Villabona
# Cross-map every subset sample's reads to every subset assembly
# ============================================================
# One array task per (assembly, reads) pair from backmap_pairs.tsv.
# This is the only step that differs in kind from the single-sample
# pipeline: instead of mapping a sample's reads only to its own contigs,
# each assembly receives reads from many samples, so the depth matrix
# built in 02 carries coverage covariance across samples.
#
# Input:  ${RESULTS}/backmap_pairs.tsv  (from 00_make_backmap_pairs.sh)
#         ${ASSEMBLY_DIR}/<asm>/final.contigs.fa          (read-only)
#         ${READS_DIR}/<rd>/<rd>_R1/R2_nonhost.fastq.gz   (read-only)
# Output: ${SCRATCH}/01_crossmap/<asm>/<asm>__<rd>.sorted.bam (+ .bai)
#         ${SCRATCH}/00_index/<asm>.mmi
#
# Tools: Minimap2 2.28 (-ax sr), Samtools 1.15.1
#        Both from modules — no conda needed in this step.
#
# Notes:
# - Pipes straight to a sorted BAM; never writes an intermediate SAM
#   (484 SAM files would be hundreds of GB of pointless I/O).
# - BAMs are indexed here because CONCOCT's coverage table step in 02
#   needs .bai files.
# - Reads are mapped WITHOUT the singleton files, matching how the
#   single-sample depth was built in 3a. Charlotte Huyghe's advice is that
#   singletons should not be used at all; they are in the assemblies
#   already, but adding them only to the multi arm would confound this.
#
# Usage: bash 00_make_backmap_pairs.sh          # writes the pair list
#        sbatch --array=1-484%20 01_crossmap.sh # N = lines in that list
# ============================================================

set -euo pipefail

# Under sbatch the script runs from SLURM's spool directory, so BASH_SOURCE does
# NOT point at this directory (first run, 2026-11-25: all 22 tasks died with
# "config.sh: No such file or directory"). Use the absolute path; override with
# SCRIPT_DIR=... if the repo ever moves.
SCRIPT_DIR="${SCRIPT_DIR:-/dfs10/hammert-lab/nvillabo/Ch2/SYMBUS/scripts/12_multisample_binning}"
source "${SCRIPT_DIR}/config.sh"

module purge
module load minimap2/2.28
module load samtools/1.15.1

THREADS="${SLURM_CPUS_PER_TASK:-8}"

: "${SLURM_ARRAY_TASK_ID:?Run with --array=1-N; N = number of lines in backmap_pairs.tsv}"
[[ -s "${PAIRS_FILE}" ]] || die "pair list not found: ${PAIRS_FILE} — run 00_make_backmap_pairs.sh first"

# ======== RESOLVE THIS TASK'S PAIR ========
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PAIRS_FILE}")
[[ -n "${LINE}" ]] || { echo "[$(date)] No pair at index ${SLURM_ARRAY_TASK_ID} — exiting."; exit 0; }
ASM=$(cut -f1 <<< "${LINE}")
RD=$(cut -f2 <<< "${LINE}")

CONTIGS="${ASSEMBLY_DIR}/${ASM}/final.contigs.fa"
R1="${READS_DIR}/${RD}/${RD}_R1_nonhost.fastq.gz"
R2="${READS_DIR}/${RD}/${RD}_R2_nonhost.fastq.gz"

IDX_DIR="${SCRATCH}/00_index"
OUT_DIR="${SCRATCH}/01_crossmap/${ASM}"
BAM="${OUT_DIR}/${ASM}__${RD}.sorted.bam"
IDX="${IDX_DIR}/${ASM}.mmi"

mkdir -p "${IDX_DIR}" "${OUT_DIR}" "${LOG_DIR}"

echo "[$(date)] task ${SLURM_ARRAY_TASK_ID}: assembly=${ASM} reads=${RD} threads=${THREADS}"

# ======== SKIP IF DONE ========
if [[ -s "${BAM}" && -s "${BAM}.bai" ]]; then
    echo "[$(date)] BAM + index already exist, skipping: ${BAM}"
    exit 0
fi

[[ -s "${CONTIGS}" ]] || die "missing assembly: ${CONTIGS}"
[[ -s "${R1}" && -s "${R2}" ]] || die "missing reads for ${RD}"

# ======== MINIMAP2 INDEX (one per assembly, shared by ~22 tasks) ========
# Tasks for the same assembly race here, so build under an O_EXCL lock and
# have the losers wait. Same pattern as 1b_host_removal.sh.
LOCK="${IDX}.lock"
if [[ ! -s "${IDX}" ]]; then
    if ( set -o noclobber; echo "$$" > "${LOCK}" ) 2>/dev/null; then
        trap 'rm -f "${LOCK}"' EXIT
        echo "[$(date)] building minimap2 index for ${ASM}"
        minimap2 -d "${IDX}.tmp.$$" "${CONTIGS}"
        mv "${IDX}.tmp.$$" "${IDX}"
        rm -f "${LOCK}"
        trap - EXIT
    else
        echo "[$(date)] another task is building ${ASM}.mmi — waiting"
        WAITED=0
        while [[ ! -s "${IDX}" ]]; do
            sleep 10
            WAITED=$((WAITED+10))
            [[ "${WAITED}" -lt 1800 ]] || die "timed out waiting for ${IDX}"
        done
    fi
fi

# ======== MAP -> SORTED BAM ========
# -T keeps samtools sort's temp files unique: several array tasks can land on
# the same node and would otherwise collide on default temp prefixes.
SORT_TMP="${OUT_DIR}/.sorttmp_${SLURM_ARRAY_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$(dirname "${SORT_TMP}")"

echo "[$(date)] mapping ${RD} -> ${ASM}"
minimap2 -ax sr --secondary=no -t "${THREADS}" "${IDX}" "${R1}" "${R2}" \
  | samtools view -@ 2 -b - \
  | samtools sort -@ "${THREADS}" -T "${SORT_TMP}" -o "${BAM}.tmp" -
mv "${BAM}.tmp" "${BAM}"
samtools index -@ "${THREADS}" "${BAM}"
rm -rf "${SORT_TMP}"*

# ======== REPORT ========
MAPPED=$(samtools view -c -F 4 "${BAM}")
TOTAL=$(samtools view -c "${BAM}")
SIZE=$(du -h "${BAM}" | cut -f1)
echo "[$(date)] DONE ${ASM}__${RD}: ${MAPPED}/${TOTAL} reads mapped, ${SIZE}"
