#!/usr/bin/env bash
#SBATCH --job-name=prodigal_snodv2
#SBATCH --output=logs/01_prodigal_%A_%a.out
#SBATCH --error=logs/01_prodigal_%A_%a.err
#SBATCH --array=0-55         # 46 MAGs + 10 refs = 56 genomes (0-indexed)
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --partition=standard

# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Gene prediction with Prodigal
# Runs on every genome in the manifest (refs + our MAGs)
# Produces .faa (protein) for OrthoFinder
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail
source ~/.bashrc

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v2"
MANIFEST="${ROOT}/00_metadata/all_genomes.txt"
mkdir -p "${ROOT}/logs" "${ROOT}/03_prodigal"

# Read genome_id and path from manifest
mapfile -t IDS  < <(cut -f1 "${MANIFEST}")
mapfile -t PATHS < <(cut -f2 "${MANIFEST}")

GENOME_ID="${IDS[$SLURM_ARRAY_TASK_ID]}"
FA="${PATHS[$SLURM_ARRAY_TASK_ID]}"
OUT="${ROOT}/03_prodigal/${GENOME_ID}"
mkdir -p "${OUT}"

echo "Genome: ${GENOME_ID}"
echo "FASTA:  ${FA}"

conda activate prodigal-env 2>/dev/null || conda activate base

prodigal \
    -i "${FA}" \
    -a "${OUT}/${GENOME_ID}.faa" \
    -d "${OUT}/${GENOME_ID}.ffn" \
    -f gff \
    -o "${OUT}/${GENOME_ID}.gff" \
    -p meta \
    -q

# Clean headers for OrthoFinder (no spaces or special chars)
sed -i 's/ /_/g; s/\t/_/g' "${OUT}/${GENOME_ID}.faa"

echo "Done: ${GENOME_ID}  ($(grep -c '^>' "${OUT}/${GENOME_ID}.faa") proteins)"
