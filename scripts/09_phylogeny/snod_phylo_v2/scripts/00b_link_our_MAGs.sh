#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# 00b_link_our_MAGs.sh
# Link our Snodgrassella MAGs into the v2 pipeline directory
# Run from LOGIN NODE
#
# Uses the 7 focal Snodgrassella MAGs from §20 of the Rmd.
# If you want to include ALL 46 MAGs from the original phylogeny,
# set USE_ALL_MAGS=true below.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v2"
MAG_SRC="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected"
OUTDIR="${ROOT}/02_our_MAGs"
mkdir -p "${OUTDIR}"

# ── Option 1: 7 focal MAGs (novel Colombian + known) ─────────────────────────
FOCAL_MAGS=(
    "NV006__54_sub.fa"
    "NV015__NV015_metabat.8.fa"
    "NV022__NV022_metabat.5.fa"
    "NV048__NV048_metabat.2.fa"
    "NV018__13_sub.fa"
    "202006191520__22_sub.fa"
    "202206071041G__39.fa"
)

# ── Option 2: ALL MAGs from original phylogeny ────────────────────────────────
# To use all 46, set USE_ALL_MAGS=true and provide the list file.
# The list file should have one MAG filename per line.
USE_ALL_MAGS=true
ALL_MAGS_LIST="${ROOT}/00_metadata/all_snod_mags.txt"  # create this if needed

# ── Copy/link ─────────────────────────────────────────────────────────────────
if [[ "${USE_ALL_MAGS}" == "true" && -f "${ALL_MAGS_LIST}" ]]; then
    echo "Using ALL MAGs from list: ${ALL_MAGS_LIST}"
    while IFS= read -r mag; do
        [[ -z "$mag" ]] && continue
        if [[ -f "${MAG_SRC}/${mag}" ]]; then
            cp "${MAG_SRC}/${mag}" "${OUTDIR}/"
            echo "  Copied: ${mag}"
        else
            echo "  MISSING: ${MAG_SRC}/${mag}"
        fi
    done < "${ALL_MAGS_LIST}"
else
    echo "Using 7 focal Snodgrassella MAGs"
    for mag in "${FOCAL_MAGS[@]}"; do
        if [[ -f "${MAG_SRC}/${mag}" ]]; then
            cp "${MAG_SRC}/${mag}" "${OUTDIR}/"
            echo "  Copied: ${mag}"
        else
            echo "  MISSING: ${MAG_SRC}/${mag}"
        fi
    done
fi

echo ""
echo "MAGs in ${OUTDIR}/:"
ls -lh "${OUTDIR}/"
