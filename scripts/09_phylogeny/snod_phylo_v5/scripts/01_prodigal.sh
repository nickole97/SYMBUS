#!/usr/bin/env bash
#SBATCH --job-name=prodigal_snodv5
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

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v5"
# v3 pointed at 06_dereplication/derep_genomes, which does not exist; the MAG
# FASTAs live in 06_MAGs_selected (140 .fa). That is why v3 silently reused
# v2's .faa copies (and inherited two Gilliamella). Fixed here.
MAG_DIR="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected"
MAG_LIST="${ROOT}/00_metadata/all_snod_mags.txt"
FAA_OUT="${ROOT}/02_prodigal"
mkdir -p "${FAA_OUT}" "${ROOT}/logs"

echo "=== Step 1: Prodigal gene prediction ==="
echo "MAG list: $(wc -l < "${MAG_LIST}") genomes"

# ── MAGs ──────────────────────────────────────────────────────────────────────
while IFS= read -r MAG; do
    FA="${MAG_DIR}/${MAG}.fa"
    if [[ ! -f "${FA}" ]]; then
        echo "ERROR: ${FA} not found (a silent skip here is how v3 lost genomes)" >&2
        exit 1
    fi
    prodigal -i "${FA}" -a "${FAA_OUT}/${MAG}.faa" -p meta -f gff -q
    echo "  MAG done: ${MAG}"
done < "${MAG_LIST}"

# ── REF + OUTGROUP genomes ────────────────────────────────────────────────────
# Reference + outgroup proteomes: reuse the Prodigal .faa produced for v3
# (same Prodigal version/settings), EXCLUDING the two Gilliamella genomes
# (GCA_029433345.1 = ESL0727, GCA_029433365.1 = ESL0721) that leaked into v3
# through .faa copied from v2. Listed explicitly so nothing else can leak.
V3_FAA="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v3/02_prodigal"
REF_FAA=(
    GCF_000600005.1   # S. alvi wkB2 (type)
    GCF_030758695.1   # S. alvi K-FP26
    GCF_002777855.1   # S. alvi, A. cerana
    GCF_002777655.1   # S. alvi, A. florea
    GCF_002777775.1   # S. alvi, A. andreniformis
    GCF_002777335.1   # S. communis type
    GCA_914768025.1   # S. gandavensis type
    GCA_029433475.1   # ESL0693, Snodgrassella from stingless bee
    Simonsiella_muelleri         # outgroup
    Rivihabitans_pingtungensis   # outgroup
)
echo "Copying ${#REF_FAA[@]} reference/outgroup proteomes from v3..."
for NAME in "${REF_FAA[@]}"; do
    [[ -s "${V3_FAA}/${NAME}.faa" ]] || { echo "ERROR: missing ${V3_FAA}/${NAME}.faa" >&2; exit 1; }
    cp "${V3_FAA}/${NAME}.faa" "${FAA_OUT}/${NAME}.faa"
    echo "  REF/OUTGROUP copied: ${NAME}"
done
N_FAA=$(ls "${FAA_OUT}"/*.faa | wc -l)
[[ "${N_FAA}" -eq $(( $(wc -l < "${MAG_LIST}") + ${#REF_FAA[@]} )) ]] || \
    { echo "ERROR: expected $(( $(wc -l < "${MAG_LIST}") + ${#REF_FAA[@]} )) .faa, found ${N_FAA}" >&2; exit 1; }
echo ""
echo "Prodigal complete."
echo "FAA files: $(ls "${FAA_OUT}"/*.faa | wc -l)"
