#!/usr/bin/env bash
#SBATCH --job-name=iqtree_snodv4
#SBATCH --output=logs/03_iqtree_%j.out
#SBATCH --error=logs/03_iqtree_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab

# ─────────────────────────────────────────────────────────────────────────────
# Step 3: IQ-TREE2 — supermatrix partitioned phylogeny.
# Uses -q (edge-unlinked) to avoid segfault with >200 partitions (IQ-TREE 2.4.0 bug).
# ─────────────────────────────────────────────────────────────────────────────

source ~/.bashrc
conda activate phylo
set -euo pipefail

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v4"
ALN="${ROOT}/02_alignments/supermatrix.faa"
PART="${ROOT}/02_alignments/partitions.txt"
OUT="${ROOT}/03_iqtree"
mkdir -p "${OUT}"

echo "Input:      ${ALN}"
echo "Partitions: ${PART}"
echo "Threads:    ${SLURM_CPUS_PER_TASK}"
echo "Taxa:       $(grep -c '^>' "${ALN}")"
echo "Partitions: $(wc -l < "${PART}")"

iqtree2 \
    -s  "${ALN}" \
    -q  "${PART}" \
    -m  TEST \
    -bb 1000 \
    -alrt 1000 \
    -nt "${SLURM_CPUS_PER_TASK}" \
    -pre "${OUT}/snod_phylo_v4" \
    --redo

echo ""
echo "IQ-TREE2 complete."
echo "Tree: ${OUT}/snod_phylo_v4.treefile"
