#!/usr/bin/env bash
#SBATCH --job-name=align_snodv4
#SBATCH --output=logs/02_align_%j.out
#SBATCH --error=logs/02_align_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab

# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Re-align augmented OG FAAs (ingroup + outgroup hits), trim, concat.
# Handles missing outgroup taxa by filling with gap characters in supermatrix.
# ─────────────────────────────────────────────────────────────────────────────

source ~/.bashrc
conda activate phylo
set -euo pipefail

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v4"
AUG_DIR="${ROOT}/01_augmented_OGs"
ALN_DIR="${ROOT}/02_alignments"
SCRIPTS="/dfs10/hammert-lab/nvillabo/Ch2/SYMBUS/scripts/09_phylogeny/snod_phylo_v4/scripts"
mkdir -p "${ALN_DIR}/raw" "${ALN_DIR}/trimmed"

N=$(ls "${AUG_DIR}"/*.faa | wc -l)
echo "Re-aligning ${N} augmented orthogroups..."

for FASTA in "${AUG_DIR}"/*.faa; do
    OG=$(basename "${FASTA}" .faa)

    mafft \
        --auto \
        --thread "${SLURM_CPUS_PER_TASK}" \
        --quiet \
        "${FASTA}" \
        > "${ALN_DIR}/raw/${OG}.aln"

    trimal \
        -in  "${ALN_DIR}/raw/${OG}.aln" \
        -out "${ALN_DIR}/trimmed/${OG}.trimmed.aln" \
        -automated1
done

echo "Alignments complete. Concatenating supermatrix..."

python3 "${SCRIPTS}/02_concat_supermatrix.py" \
    --trimmed_dir "${ALN_DIR}/trimmed" \
    --outdir      "${ALN_DIR}"

echo ""
echo "Supermatrix stats:"
grep -c "^>" "${ALN_DIR}/supermatrix.faa" | xargs -I{} echo "  Taxa: {}"
awk '/^>/{next}{len+=length($0)}END{print "  Length:", len, "aa"}' \
    <(awk '/^>/{if(seq)print seq; seq=""; print; next}{seq=seq $0}END{if(seq)print seq}' \
      "${ALN_DIR}/supermatrix.faa") 2>/dev/null || true
echo "Partition file: ${ALN_DIR}/partitions.txt"
