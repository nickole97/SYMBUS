#!/usr/bin/env bash
#SBATCH --job-name=host_rm
#SBATCH --output=host_rm.o%A_%a
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G

# ============================================================
# 1b_host_removal.sh
# Author: Nickole Villabona
# Remove host (Bombus impatiens) reads with Minimap2
# ============================================================
# Maps cleaned reads against the bee host genome and filters
# out mapped reads, keeping only non-host (microbial) reads.
#
# Input:  /work_meta/preprocessing/01_preprocessing/<sample>/
#         - *_R1_clean.fastq.gz, *_R2_clean.fastq.gz
# Output: /work_meta/preprocessing/02_host_removal/<sample>/
#         - *_R1_nonhost.fastq.gz, *_R2_nonhost.fastq.gz
#         - *_singleton_nonhost.fastq.gz
#         - *_combined_nonhost.fastq.gz (all three concatenated)
#
# Reference: Bombus impatiens BIMP_2.2 (GCA_000188095.4)
# Tools: Minimap2 (mapping), Samtools (filtering unmapped reads)
#
# Usage: sbatch --array=0-112 1b_host_removal.sh
# ============================================================

set -euo pipefail

# ====== CONFIG ======
BASE_WD="/dfs10/hammert-lab/nvillabo/Ch2/work_meta"
WD="${BASE_WD}/preprocessing"

# Host genome (Bombus impatiens BIMP_2.2)
REF="${BASE_WD}/GCA_000188095.4_BIMP_2.2_genomic.fna"
REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/188/095/GCA_000188095.4_BIMP_2.2/GCA_000188095.4_BIMP_2.2_genomic.fna.gz"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_NUM_THREADS="${THREADS}"

IDXDIR="${WD}/host_index"
IDX="${IDXDIR}/host.mmi"
SAMPLES="${WD}/samples.txt"

# ====== MODULES ======
module purge
module load minimap2/2.28
module load samtools/1.15.1
module load pigz/2.6 || true

# ====== PREP ======
mkdir -p "${WD}/"{logs,02_host_removal} "${IDXDIR}"

# Download host reference if missing
if [[ ! -f "${REF}" ]]; then
  echo "[REF] Downloading reference…"
  mkdir -p "${BASE_WD}" && cd "${BASE_WD}"
  wget -q "${REF_URL}" -O hostref.fna.gz
  gunzip -f hostref.fna.gz && mv hostref.fna "${REF}"
fi

# Minimap2 index (simple lock to avoid race conditions)
LOCK="${IDX}.lock"
if [[ ! -f "${IDX}" ]]; then
  if ( set -o noclobber; echo "$$" > "${LOCK}" ) 2>/dev/null; then
    trap 'rm -f "${LOCK}"' EXIT
    echo "[IDX] Building index…"
    minimap2 -d "${IDX}" "${REF}" 2>&1 | tee "${WD}/logs/index_host.log"
  else
    echo "[IDX] Another task is building the index; waiting…"
    while [[ ! -f "${IDX}" ]]; do sleep 5; done
  fi
fi

# Sample list (generated once from clean reads)
if [[ ! -s "${SAMPLES}" ]]; then
  echo "[SAMPLES] Generating list…"
  find "${WD}/01_preprocessing" -mindepth 1 -maxdepth 1 -type d \
    -exec bash -c 'd="$1"; sid=$(basename "$d"); [[ -s "$d/${sid}_R1_clean.fastq.gz" && -s "$d/${sid}_R2_clean.fastq.gz" ]] && echo "$sid"' _ {} \; \
    | sort > "${SAMPLES}"
fi

# ====== ARRAY SELECTION ======
: "${SLURM_ARRAY_TASK_ID:?Run with --array=0-(N-1)}"
SID="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${SAMPLES}" || true)"
[[ -n "${SID}" ]] || { echo "[ARRAY] No sample for index ${SLURM_ARRAY_TASK_ID}"; exit 0; }

PRE="${WD}/01_preprocessing/${SID}"
OUT="${WD}/02_host_removal/${SID}"
LOG="${WD}/logs/host_${SID}.log"
mkdir -p "${OUT}"

echo "[TASK] ${SID} | thr=${THREADS} | $(date)"

# Skip if combined file already exists
if [[ -s "${OUT}/${SID}_combined_nonhost.fastq.gz" ]]; then
  echo "✅ Combined file already exists: ${SID}"; exit 0
fi

# Check inputs
[[ -s "${PRE}/${SID}_R1_clean.fastq.gz" && -s "${PRE}/${SID}_R2_clean.fastq.gz" ]] || { echo "⚠️ Missing clean R1/R2"; exit 0; }

# Parallel compressor
if command -v pigz >/dev/null 2>&1; then
  COMP=(pigz -p "${THREADS}" -1)
else
  COMP=(gzip -1)
fi

# ====== MAPPING + FILTERING + FASTQ ======
set -o pipefail
minimap2 -ax sr --secondary=no -K 2g -t "${THREADS}" "${IDX}" \
  "${PRE}/${SID}_R1_clean.fastq.gz" "${PRE}/${SID}_R2_clean.fastq.gz" 2> "${LOG}" \
| samtools view -@ "${THREADS}" -b -F 256 -f 4 - \
| samtools fastq -@ "${THREADS}" \
    -1 >("${COMP[@]}" > "${OUT}/${SID}_R1_nonhost.fastq.gz") \
    -2 >("${COMP[@]}" > "${OUT}/${SID}_R2_nonhost.fastq.gz") \
    -s >("${COMP[@]}" > "${OUT}/${SID}_singleton_nonhost.fastq.gz") \
    - 2>> "${LOG}"

# ====== COMBINE ======
cat "${OUT}/${SID}_R1_nonhost.fastq.gz" \
    "${OUT}/${SID}_R2_nonhost.fastq.gz" \
    "${OUT}/${SID}_singleton_nonhost.fastq.gz" \
  > "${OUT}/${SID}_combined_nonhost.fastq.gz"

echo "[DONE] ${SID} | $(date)"
