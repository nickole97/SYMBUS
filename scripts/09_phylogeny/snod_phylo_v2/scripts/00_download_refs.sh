#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# 00_download_refs.sh
# Download Snodgrassella reference genomes from NCBI using datasets CLI
# Run from LOGIN NODE — only network transfer, no computation
#
# BEFORE RUNNING: confirm GCA accessions for ESL0693/ESL0721/ESL0727 at:
#   https://www.ncbi.nlm.nih.gov/assembly/?term=ESL0693
#   https://www.ncbi.nlm.nih.gov/assembly/?term=ESL0721
#   https://www.ncbi.nlm.nih.gov/assembly/?term=ESL0727
# Then fill in STINGLESS_ACCS below.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

ROOT="$HOME/snod_phylo_v2"
OUTDIR="${ROOT}/01_reference_genomes"
mkdir -p "${OUTDIR}"

# ── Check that ncbi datasets CLI is available ─────────────────────────────────
if ! command -v datasets &>/dev/null; then
    echo "ERROR: 'datasets' CLI not found."
    echo "Install with: conda install -c conda-forge ncbi-datasets-cli"
    exit 1
fi

# ── Accessions to download ────────────────────────────────────────────────────

# S. alvi from Apis (confirmed accessions)
APIS_ACCS=(
    "GCF_000600005.1"   # S. alvi wkB2 — type strain — A. mellifera USA
    "CP132376"          # S. alvi K-FP26 — A. mellifera iberiensis Portugal
    "GCF_002777855.1"   # novel Apis gr.1 — A. cerana Singapore
    "GCF_002777655.1"   # novel Apis gr.1 — A. florea Singapore
    "GCF_002777775.1"   # novel Apis gr.1 — A. andreniformis Singapore
)

# Named Bombus species type strains
BOMBUS_TYPE_ACCS=(
    "GCF_002777335.1"   # S. communis type — B. griseocollis USA
    "GCA_914768025.1"   # S. gandavensis type — B. pascuorum Belgium
)

# Stingless bee strains — FILL IN AFTER CHECKING NCBI
# TaxID 2983239 → ESL0693 (Neisseriaceae/Snodgrassella)
# TaxID 2983241 → ESL0721 (Orbaceae/Gilliamella)
# TaxID 2983240 → ESL0727 (Orbaceae/Gilliamella)
STINGLESS_ACCS=(
    "GCA_029433475.1"   # ESL0693 — Snodgrassella from stingless bee — TaxID 2983239
    "GCA_029433365.1"   # ESL0721 — Gilliamella from stingless bee — TaxID 2983241
    "GCA_029433345.1"   # ESL0727 — Gilliamella from stingless bee — TaxID 2983240
)

ALL_ACCS=("${APIS_ACCS[@]}" "${BOMBUS_TYPE_ACCS[@]}" "${STINGLESS_ACCS[@]}")

# ── Download ──────────────────────────────────────────────────────────────────
echo "Downloading ${#ALL_ACCS[@]} reference genomes..."

datasets download genome accession \
    "${ALL_ACCS[@]}" \
    --include genome \
    --filename "${OUTDIR}/refs_download.zip"

# Unzip and reorganize
cd "${OUTDIR}"
unzip -o refs_download.zip -d refs_unpacked/

# Move each .fna to a clean filename: {accession}.fna
find refs_unpacked/ -name "*.fna" | while read f; do
    acc=$(basename "$(dirname "$f")")
    cp "$f" "${OUTDIR}/${acc}.fna"
    echo "  Saved: ${acc}.fna"
done

echo ""
echo "Done. Reference genomes in: ${OUTDIR}/"
ls -lh "${OUTDIR}"/*.fna 2>/dev/null || echo "(no .fna files yet — check download)"
