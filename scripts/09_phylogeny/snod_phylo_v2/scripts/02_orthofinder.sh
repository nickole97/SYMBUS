#!/usr/bin/env bash
#SBATCH --job-name=orthofinder_snodv2
#SBATCH --output=logs/02_orthofinder_%j.out
#SBATCH --error=logs/02_orthofinder_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --partition=standard

# ─────────────────────────────────────────────────────────────────────────────
# Step 2: OrthoFinder — orthogroups across all Snodgrassella genomes
# Collects all .faa from 03_prodigal/ and runs OrthoFinder
# ─────────────────────────────────────────────────────────────────────────────

source ~/.bashrc
conda activate orthofinder-env 2>/dev/null || conda activate base
set -euo pipefail

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v2"
MANIFEST="${ROOT}/00_metadata/all_genomes.txt"
FAA_INPUT="${ROOT}/04_orthofinder/faa_input"
mkdir -p "${FAA_INPUT}" "${ROOT}/logs"

# Collect .faa files into a single input directory
while IFS=$'\t' read -r genome_id path; do
    faa="${ROOT}/03_prodigal/${genome_id}/${genome_id}.faa"
    if [[ -f "${faa}" ]]; then
        cp "${faa}" "${FAA_INPUT}/${genome_id}.faa"
    else
        echo "WARNING: missing .faa for ${genome_id} at ${faa}"
    fi
done < "${MANIFEST}"

N=$(ls "${FAA_INPUT}"/*.faa | wc -l)
echo "Input FAAs: ${N}"
ls -lh "${FAA_INPUT}/"

# Run OrthoFinder
orthofinder \
    -f "${FAA_INPUT}" \
    -t ${SLURM_CPUS_PER_TASK} \
    -a ${SLURM_CPUS_PER_TASK} \
    -o "${ROOT}/04_orthofinder/results" \
    -M msa \
    -ot

echo ""
echo "OrthoFinder complete."
echo "Results in: ${ROOT}/04_orthofinder/results"
