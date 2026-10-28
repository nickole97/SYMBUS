#!/usr/bin/env bash
#SBATCH --job-name=iqtree_snodv2
#SBATCH --output=logs/05_iqtree_%j.out
#SBATCH --error=logs/05_iqtree_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab

# ─────────────────────────────────────────────────────────────────────────────
# Step 5: IQ-TREE2 — maximum-likelihood phylogeny from protein supermatrix
# -m TEST    → model selection per partition
# -bb 1000   → ultrafast bootstrap
# -alrt 1000 → SH-aLRT branch support
# ─────────────────────────────────────────────────────────────────────────────

source ~/.bashrc
conda activate phylo
set -euo pipefail

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v2"
ALN="${ROOT}/06_alignments/supermatrix.faa"
PART="${ROOT}/06_alignments/partitions.txt"
OUT="${ROOT}/07_iqtree"
mkdir -p "${OUT}" "${ROOT}/logs"

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
    -pre "${OUT}/snod_phylo_v2" \
    --redo

echo ""
echo "IQ-TREE2 complete."
echo "Tree: ${OUT}/snod_phylo_v2.treefile"
