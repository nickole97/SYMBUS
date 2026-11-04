#!/usr/bin/env bash
#SBATCH --job-name=prodigal_gilliv2
#SBATCH --output=logs/01_prodigal_%j.out
#SBATCH --error=logs/01_prodigal_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab

source ~/.bashrc
conda activate gtdbtk-2.5.2
set -euo pipefail

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/gilli_phylo_v2"
MAG_DIR="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected"
MAG_LIST="${ROOT}/00_metadata/all_gilli_mags.txt"
REF_DIR="${ROOT}/01_genomes/refs"
FAA_OUT="${ROOT}/02_prodigal"
mkdir -p "${FAA_OUT}" "${ROOT}/logs"

echo "=== Step 1: Prodigal gene prediction ==="
echo "MAG list: $(wc -l < "${MAG_LIST}") genomes"

# ── MAGs ──────────────────────────────────────────────────────────────────────
while IFS= read -r MAG; do
    FA="${MAG_DIR}/${MAG}.fa"
    if [[ ! -f "${FA}" ]]; then
        echo "WARNING: ${FA} not found — skipping" >&2
        continue
    fi
    prodigal -i "${FA}" -a "${FAA_OUT}/${MAG}.faa" -p meta -f gff -q
    echo "  MAG done: ${MAG}"
done < "${MAG_LIST}"

# ── REF + OUTGROUP genomes ────────────────────────────────────────────────────
echo "Running Prodigal on reference/outgroup genomes..."
for FNA in "${REF_DIR}"/*.fna; do
    NAME=$(basename "${FNA}" .fna)
    prodigal -i "${FNA}" -a "${FAA_OUT}/${NAME}.faa" -p single -f gff -q
    echo "  REF/OUTGROUP done: ${NAME}"
done

echo ""
echo "Prodigal complete."
echo "FAA files: $(ls "${FAA_OUT}"/*.faa | wc -l)"
