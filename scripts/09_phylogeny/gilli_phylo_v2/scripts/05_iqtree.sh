#!/usr/bin/env bash
#SBATCH --job-name=iqtree_gilliv2
#SBATCH --output=logs/05_iqtree_%j.out
#SBATCH --error=logs/05_iqtree_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab

source ~/.bashrc
conda activate phylo
set -euo pipefail

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/gilli_phylo_v2"
ALN="${ROOT}/05_alignments/supermatrix.faa"
PART="${ROOT}/05_alignments/partitions.txt"
OUT="${ROOT}/06_iqtree"
mkdir -p "${OUT}"

echo "Input:      ${ALN}"
echo "Partitions: ${PART}"
echo "Threads:    ${SLURM_CPUS_PER_TASK}"

iqtree2 \
    -s "${ALN}" \
    -q "${PART}" \
    -m TEST \
    -bb 1000 \
    -alrt 1000 \
    -nt ${SLURM_CPUS_PER_TASK} \
    -pre "${OUT}/gilli_phylo_v2" \
    --redo

echo ""
echo "IQ-TREE2 complete."
echo "Tree: ${OUT}/gilli_phylo_v2.treefile"
