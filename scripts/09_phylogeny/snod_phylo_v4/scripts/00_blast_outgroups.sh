#!/usr/bin/env bash
#SBATCH --job-name=blast_outgroups_snodv4
#SBATCH --output=logs/00_blast_%j.out
#SBATCH --error=logs/00_blast_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab

# ─────────────────────────────────────────────────────────────────────────────
# Step 0: BLAST outgroup FAAs against each core OG to find orthologs.
# Uses ingroup-only core OGs from snod_phylo_v2 (411 genes).
# For each OG × outgroup: returns best-hit protein ID and sequence.
# ─────────────────────────────────────────────────────────────────────────────

source ~/.bashrc
conda activate phylo
set -euo pipefail

V2="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v2"
V3="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v3"
ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v4"

OG_DIR="${V2}/05_core_OGs"
OG_LIST="${OG_DIR}/core_OG_ids.txt"
OUT_FAA_DIR="${V3}/02_prodigal"       # Simonsiella_muelleri.faa + Rivihabitans_pingtungensis.faa
BLAST_DIR="${ROOT}/00_blast"
mkdir -p "${BLAST_DIR}" "${ROOT}/logs"

OUTGROUPS=("Simonsiella_muelleri" "Rivihabitans_pingtungensis")
EVALUE="1e-5"
THREADS="${SLURM_CPUS_PER_TASK}"

echo "Core OGs: $(wc -l < "${OG_LIST}")"
echo "Outgroups: ${OUTGROUPS[*]}"

for OG_NAME in "${OUTGROUPS[@]}"; do
    FAA="${OUT_FAA_DIR}/${OG_NAME}.faa"
    DB_DIR="${BLAST_DIR}/${OG_NAME}_db"
    mkdir -p "${DB_DIR}"

    echo ""
    echo "=== Building BLAST DB for ${OG_NAME} ==="
    makeblastdb -in "${FAA}" -dbtype prot -out "${DB_DIR}/${OG_NAME}" -logfile /dev/null

    HIT_DIR="${BLAST_DIR}/${OG_NAME}_hits"
    mkdir -p "${HIT_DIR}"

    echo "BLASTing ${OG_NAME} against $(wc -l < "${OG_LIST}") core OGs..."
    while IFS= read -r OG; do
        QUERY="${OG_DIR}/${OG}.faa"
        OUT_TSV="${HIT_DIR}/${OG}.tsv"
        [[ -f "${OUT_TSV}" ]] && continue   # skip if already done

        blastp \
            -query     "${QUERY}" \
            -db        "${DB_DIR}/${OG_NAME}" \
            -out       "${OUT_TSV}" \
            -outfmt    "6 qseqid sseqid pident length evalue bitscore" \
            -evalue    "${EVALUE}" \
            -max_hsps  1 \
            -max_target_seqs 1 \
            -num_threads "${THREADS}" \
            2>/dev/null
    done < "${OG_LIST}"

    NHITS=$(find "${HIT_DIR}" -name "*.tsv" -size +0 | wc -l)
    echo "  OGs with outgroup hit: ${NHITS} / $(wc -l < "${OG_LIST}")"
done

echo ""
echo "BLAST complete. Results in: ${BLAST_DIR}"
