#!/usr/bin/env bash
# Download reference + outgroup genomes for gilli_phylo_v2
# Run from login node (no compute needed)
set -euo pipefail

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/gilli_phylo_v2"
mkdir -p "${ROOT}/01_genomes/refs"
cd "${ROOT}/01_genomes/refs"

declare -A GENOMES=(
    # REF ingroup — Bombus Gilliamella isolates
    [GCA_001690335.1]="Gris1-4"
    [GCA_001690595.1]="Choc5-1"
    [GCA_001690655.1]="Occ3-1"
    [GCA_013344985.1]="ESL0254"
    [GCA_013345045.1]="ESL0250"
    [GCA_013346885.1]="ESL0232"
    [GCA_900094935.1]="R-53144"
    [GCA_900094945.1]="R-53248"
    [GCA_029433365.1]="ESL0721"
    [GCA_029433345.1]="ESL0727"
    # OUTGROUP — Bombus pratorum Gilliamella (near outgroup)
    [GCA_009795805.1]="Gilli_Pra-s52"
    [GCA_009795785.1]="Gilli_Pra-s54"
    [GCA_009795755.1]="Gilli_Pra-s60"
    [GCA_009795745.1]="Gilli_Pra-s65"
    # OUTGROUP — deep root (Pasteurellaceae)
    [GCA_900186835.1]="Pasteurella_dagmatis"
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
