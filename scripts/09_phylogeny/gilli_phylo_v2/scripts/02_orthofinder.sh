#!/usr/bin/env bash
#SBATCH --job-name=orthofinder_gilliv2
#SBATCH --output=logs/02_orthofinder_%j.out
#SBATCH --error=logs/02_orthofinder_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab

source ~/.bashrc
conda activate phylo
set -euo pipefail

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/gilli_phylo_v2"
FAA_DIR="${ROOT}/02_prodigal"
OF_INPUT="${ROOT}/03_orthofinder/faa_input"
OF_OUT="${ROOT}/03_orthofinder/results"
mkdir -p "${OF_INPUT}" "${ROOT}/logs"

# Copy all FAAs into OrthoFinder input dir (flat, one .faa per genome)
cp "${FAA_DIR}"/*.faa "${OF_INPUT}/"

N=$(ls "${OF_INPUT}"/*.faa | wc -l)
echo "Input FAAs: ${N} genomes"
ls "${OF_INPUT}/"

orthofinder \
    -f "${OF_INPUT}" \
    -t ${SLURM_CPUS_PER_TASK} \
    -a ${SLURM_CPUS_PER_TASK} \
    -o "${OF_OUT}" \
    -M msa

echo ""
echo "OrthoFinder complete."
echo "Results: ${OF_OUT}"
