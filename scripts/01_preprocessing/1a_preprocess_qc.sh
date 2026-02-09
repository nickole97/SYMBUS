pero#!/usr/bin/env bash
#SBATCH --job-name=qual_pre_only
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs/qual_pre_only_%j.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs/qual_pre_only_%j.err
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

# ============================================================
# 1a_preprocess_qc.sh
# Author: Nickole Villabona
# Quality control and adapter trimming with BBDuk + BBMerge
# ============================================================
# Takes raw paired reads, removes adapters and low-quality bases,
# then merges overlapping pairs where possible.
#
# Input:  /Ch2/raw_01/raw_paired/<sample>/<sample>_{1,2}.fq.gz
# Output: /work_meta/preprocessing/01_preprocessing/<sample>/
#         - *_R1_clean.fastq.gz, *_R2_clean.fastq.gz
#         - *_merged.fastq.gz, *_R1/R2_unmerged.fastq.gz
#
# Tools: BBDuk (ktrim=r, k=23, trimq=20, minlen=50)
#        BBMerge (default parameters)
# ============================================================

set -euo pipefail

# ======== CONFIG ========
READS_DIR="/dfs10/hammert-lab/nvillabo/Ch2/raw_01/raw_paired"   # SAMPLE/SAMPLE_{1,2}.fq.gz
BASE_WD="/dfs10/hammert-lab/nvillabo/Ch2/work_meta"
WD="${BASE_WD}/preprocessing"
THREADS_PRE=4

# ======== MODULES ========
module purge
module load bbmap/38.96

echo "[INFO] QC only (bbduk+bbmerge) | READS_DIR=$READS_DIR | THREADS=$THREADS_PRE | $(date)"

# ======== PREP ========
mkdir -p "${WD}"/{logs,01_preprocessing}

# ======== SAMPLE LIST ========
mapfile -t sample_ids < <(find "$READS_DIR" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort)
echo "[IDX] ${#sample_ids[@]} samples: ${sample_ids[*]}"

# ======== FUNCTIONS ========
preprocess_reads() {
  local sid="$1"
  local r1 r2 outdir
  r1=$(ls "${READS_DIR}/${sid}/${sid}"*_1.f*q.gz 2>/dev/null | head -1 || true)
  r2=$(ls "${READS_DIR}/${sid}/${sid}"*_2.f*q.gz 2>/dev/null | head -1 || true)
  outdir="${WD}/01_preprocessing/${sid}"

  if [[ -z "${r1}" || -z "${r2}" ]]; then
    echo "[WARN] Missing R1/R2 for ${sid}"; return 1
  fi
  if [[ -s "${outdir}/${sid}_R1_clean.fastq.gz" && -s "${outdir}/${sid}_R2_clean.fastq.gz" ]]; then
    echo "✅ Preprocessing already done: ${sid}"; return 0
  fi

  mkdir -p "$outdir"
  echo "[BBDUK] ${sid}"
  bbduk.sh \
    in1="${r1}" in2="${r2}" \
    out1="${outdir}/${sid}_R1_clean.fastq.gz" \
    out2="${outdir}/${sid}_R2_clean.fastq.gz" \
    ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20 minlen=50 \
    threads=${THREADS_PRE} 2>&1 | tee "${WD}/logs/pre_${sid}.log"

  echo "[BBMERGE] ${sid}"
  bbmerge.sh \
    in1="${outdir}/${sid}_R1_clean.fastq.gz" \
    in2="${outdir}/${sid}_R2_clean.fastq.gz" \
    out="${outdir}/${sid}_merged.fastq.gz" \
    outu1="${outdir}/${sid}_R1_unmerged.fastq.gz" \
    outu2="${outdir}/${sid}_R2_unmerged.fastq.gz" \
    threads=${THREADS_PRE} 2>&1 | tee -a "${WD}/logs/pre_${sid}.log"
}

# ======== EXECUTION ========
for sid in "${sample_ids[@]}"; do
  preprocess_reads "$sid"
done

echo "[DONE] QC finished. Clean reads in: ${WD}/01_preprocessing/<SAMPLE>/  | $(date)"
