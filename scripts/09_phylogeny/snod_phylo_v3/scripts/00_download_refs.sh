#!/usr/bin/env bash
# Download reference + outgroup genomes for snod_phylo_v3
# Run from login node (no compute needed — just wget/datasets)
set -euo pipefail

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v3"
mkdir -p "${ROOT}/01_genomes/refs"
cd "${ROOT}/01_genomes/refs"

# ── Tool: NCBI datasets CLI (available on HPC3) ───────────────────────────────
# If not available: module load ncbi-datasets or conda activate phylo

declare -A GENOMES=(
    # REF ingroup
    [GCF_000600005.1]="S_alvi_wkB2"
    [GCF_030758695.1]="S_alvi_K-FP26"
    [GCF_002777855.1]="S_alvi_Acerana"
    [GCF_002777655.1]="S_alvi_Aflorea"
    [GCF_002777775.1]="S_alvi_Aandreniformis"
    [GCF_002777335.1]="S_communis_type"
    [GCA_914768025.1]="S_gandavensis_type"
    [GCA_029433475.1]="ESL0693"
    # OUTGROUP
    [GCA_000163775.2]="Simonsiella_muelleri"
    [GCA_003201855.1]="Rivihabitans_pingtungensis"
)

for ACC in "${!GENOMES[@]}"; do
    NAME="${GENOMES[$ACC]}"
    echo "Downloading ${ACC} (${NAME})..."
    datasets download genome accession "${ACC}" \
        --include genome \
        --filename "${ACC}.zip"
    unzip -o "${ACC}.zip" -d "${ACC}_tmp"
    FA=$(find "${ACC}_tmp" -name "*.fna" | head -1)
    cp "${FA}" "${NAME}.fna"
    rm -rf "${ACC}.zip" "${ACC}_tmp"
    echo "  → ${NAME}.fna"
done

echo ""
echo "All reference/outgroup genomes downloaded to: ${ROOT}/01_genomes/refs/"
ls -lh *.fna
