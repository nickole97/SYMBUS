#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# 00c_make_manifest.sh
# Build a manifest of all genome FASTA files (refs + our MAGs) for the pipeline.
# Run on LOGIN NODE after 00_download_refs.sh and 00b_link_our_MAGs.sh.
#
# Output: 00_metadata/all_genomes.txt
#   Format: genome_id <TAB> fasta_path
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v2"
MANIFEST="${ROOT}/00_metadata/all_genomes.txt"

: > "${MANIFEST}"

# Reference genomes (.fna)
for f in "${ROOT}/01_reference_genomes/"*.fna; do
    [[ -f "$f" ]] || continue
    genome_id=$(basename "$f" .fna)
    echo -e "${genome_id}\t${f}" >> "${MANIFEST}"
done

# Our MAGs (.fa or .fna)
for f in "${ROOT}/02_our_MAGs/"*.fa "${ROOT}/02_our_MAGs/"*.fna; do
    [[ -f "$f" ]] || continue
    ext="${f##*.}"
    genome_id=$(basename "$f" ."${ext}")
    echo -e "${genome_id}\t${f}" >> "${MANIFEST}"
done

N=$(wc -l < "${MANIFEST}")
echo "Manifest written: ${MANIFEST}"
echo "Total genomes: ${N}"
echo ""
cat "${MANIFEST}"
echo ""
echo "Use array 0-$((N-1)) in 01_prodigal.sh"
echo "Update: #SBATCH --array=0-$((N-1))"
