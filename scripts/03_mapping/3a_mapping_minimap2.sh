#!/bin/bash
#SBATCH --account=hammert_lab
#SBATCH --job-name=map_array
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/map_array_%A_%a.out
#SBATCH --partition=standard

# ============================================================
# 3a_mapping_minimap2.sh
# Author: Nickole Villabona
# Map non-host reads back to assembled contigs
# ============================================================
# Maps cleaned non-host reads back to the assembled contigs
# for each sample to calculate contig coverage and depth.
# Creates sorted BAM files and depth tables for binning.
#
# Input:  /work_meta/preprocessing/02_host_removal/<sample>/
#         - *_R1_nonhost.fastq.gz, *_R2_nonhost.fastq.gz
#         /work_meta/03_assembly/<sample>/
#         - final.contigs.fa
# Output: /work_meta/04_mapping/<sample>/
#         - <sample>.sorted.bam (+ index)
#         - <sample>_depth.txt
#
# Tools: Minimap2 (short read mapping mode)
#        Samtools (BAM sorting and indexing)
#        jgi_summarize_bam_contig_depths (from MetaBAT2)
#
# Usage: sbatch --array=1-113%5 3a_mapping_minimap2.sh
# ============================================================

set -euo pipefail

WORK_BASE=/dfs10/hammert-lab/nvillabo/Ch2/work_meta
READS_DIR=${WORK_BASE}/preprocessing/02_host_removal
ASSEMBLY_DIR=${WORK_BASE}/03_assembly
MAP_DIR=${WORK_BASE}/04_mapping

mkdir -p "${MAP_DIR}" "${WORK_BASE}/logs"

# ======== GET SAMPLE FROM ARRAY INDEX ========
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${WORK_BASE}/preprocessing/samples.txt")

if [ -z "${SAMPLE}" ]; then
    echo "[$(date)] No sample for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}, exiting."
    exit 0
fi

echo "[$(date)] Starting mapping for ${SAMPLE}"

READ1=${READS_DIR}/${SAMPLE}/${SAMPLE}_R1_nonhost.fastq.gz
READ2=${READS_DIR}/${SAMPLE}/${SAMPLE}_R2_nonhost.fastq.gz
CONTIGS=${ASSEMBLY_DIR}/${SAMPLE}/final.contigs.fa

if [ ! -f "${READ1}" ] || [ ! -f "${READ2}" ]; then
    echo "[$(date)] Missing reads for ${SAMPLE}. Skipping."
    exit 0
fi

if [ ! -f "${CONTIGS}" ]; then
    echo "[$(date)] No assembly for ${SAMPLE}. Skipping."
    exit 0
fi

# ======== MODULES ========
module load miniconda3/23.5.2
conda activate quast_env

module load minimap2/2.28
module load samtools/1.15.1

MAP_SAMPLE_DIR=${MAP_DIR}/${SAMPLE}
mkdir -p "${MAP_SAMPLE_DIR}"
cd "${MAP_SAMPLE_DIR}"

# ======== MINIMAP2 INDEX ========
if [ ! -f "${SAMPLE}.mmi" ]; then
    echo "[$(date)] Creating minimap2 index..."
    minimap2 -d "${SAMPLE}.mmi" "${CONTIGS}"
else
    echo "[$(date)] Index exists, skipping."
fi

# ======== MAPPING ========
if [ ! -f "${SAMPLE}.sorted.bam" ]; then
    echo "[$(date)] Mapping reads..."
    minimap2 -ax sr \
      "${SAMPLE}.mmi" \
      "${READ1}" \
      "${READ2}" \
      -t "${SLURM_CPUS_PER_TASK}" \
      > "${SAMPLE}.sam"

    echo "[$(date)] Converting SAM → BAM..."
    samtools view -b "${SAMPLE}.sam" | samtools sort -o "${SAMPLE}.sorted.bam"
    samtools index "${SAMPLE}.sorted.bam"
    rm "${SAMPLE}.sam"
else
    echo "[$(date)] sorted.bam already exists, skipping mapping."
fi

# ======== DEPTH CALCULATION ========
DEPTH_FILE="${MAP_SAMPLE_DIR}/${SAMPLE}_depth.txt"

if [ ! -f "${DEPTH_FILE}" ]; then
    echo "[$(date)] Calculating depth..."
    jgi_summarize_bam_contig_depths \
      --outputDepth "${DEPTH_FILE}" \
      "${SAMPLE}.sorted.bam"
else
    echo "[$(date)] Depth file already exists, skipping."
fi

echo "[$(date)] DONE mapping for ${SAMPLE}"
