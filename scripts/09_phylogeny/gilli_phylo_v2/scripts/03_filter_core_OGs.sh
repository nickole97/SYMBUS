#!/usr/bin/env bash
#SBATCH --job-name=filter_OGs_gilliv2
#SBATCH --output=logs/03_filter_%j.out
#SBATCH --error=logs/03_filter_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab

source ~/.bashrc
conda activate phylo
set -euo pipefail

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/gilli_phylo_v2"
OF_RES="${ROOT}/03_orthofinder/results"

# OrthoFinder writes results to a date-stamped subfolder — find the latest one
OF_DIR=$(ls -td "${OF_RES}"/Results_* 2>/dev/null | head -1)
if [[ -z "${OF_DIR}" ]]; then
    echo "ERROR: no OrthoFinder results found in ${OF_RES}" >&2
    exit 1
fi
echo "OrthoFinder results: ${OF_DIR}"

OG_TSV="${OF_DIR}/Orthogroups/Orthogroups.tsv"
FAA_DIR="${OF_DIR}/Orthogroup_Sequences"
OUT="${ROOT}/04_core_OGs"
mkdir -p "${OUT}"

SCRIPTS="/dfs10/hammert-lab/nvillabo/Ch2/SYMBUS/scripts/09_phylogeny/gilli_phylo_v2/scripts"
python3 "${SCRIPTS}/03_filter_core_OGs.py" \
    --orthogroups "${OG_TSV}" \
    --faa_dir     "${FAA_DIR}" \
    --outdir      "${OUT}"

echo ""
echo "Core OGs: $(wc -l < "${OUT}/core_OG_ids.txt")"
echo "Output: ${OUT}"
