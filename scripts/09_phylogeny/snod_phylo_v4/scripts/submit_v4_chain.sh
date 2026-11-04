#!/usr/bin/env bash
# Submit snod_phylo_v4 job chain with SLURM dependencies.
# Run from: /dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v4/
# Prerequisite: v2 core OGs exist at ../snod_phylo_v2/05_core_OGs/
#               v3 outgroup FAAs exist at ../snod_phylo_v3/02_prodigal/

set -euo pipefail

SCRIPTS="/dfs10/hammert-lab/nvillabo/Ch2/SYMBUS/scripts/09_phylogeny/snod_phylo_v4/scripts"
mkdir -p logs

echo "=== Submitting snod_phylo_v4 job chain ==="

# Job 0: BLAST outgroups against all 411 core OGs
JID0=$(sbatch --parsable "${SCRIPTS}/00_blast_outgroups.sh")
echo "  Job 0 (BLAST outgroups):       ${JID0}"

# Job 1: Parse BLAST results, add outgroup seqs to OG FAAs
JID1=$(sbatch --parsable \
    --dependency=afterok:${JID0} \
    --wrap="
source ~/.bashrc && conda activate phylo && set -euo pipefail

V2='/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v2'
ROOT='/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v4'
SCRIPTS='/dfs10/hammert-lab/nvillabo/Ch2/SYMBUS/scripts/09_phylogeny/snod_phylo_v4/scripts'

mkdir -p \"\${ROOT}/01_augmented_OGs\" \"\${ROOT}/logs\"

python3 \"\${SCRIPTS}/01_add_outgroup_seqs.py\" \
    --v2_og_dir   \"\${V2}/05_core_OGs\" \
    --blast_dir   \"\${ROOT}/00_blast\" \
    --out_faa_dir '/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v3/02_prodigal' \
    --outdir      \"\${ROOT}/01_augmented_OGs\" \
    --evalue      1e-5
" \
    --job-name=add_outgroup_snodv4 \
    --output=logs/01_add_outgroup_%j.out \
    --error=logs/01_add_outgroup_%j.err \
    --cpus-per-task=2 \
    --mem=8G \
    --time=01:00:00 \
    --partition=standard \
    --account=hammert_lab)
echo "  Job 1 (add outgroup seqs):     ${JID1}"

# Job 2: Re-align with MAFFT, trim with trimAl, concat supermatrix
JID2=$(sbatch --parsable --dependency=afterok:${JID1} "${SCRIPTS}/02_realign_concat.sh")
echo "  Job 2 (align + concat):        ${JID2}"

# Job 3: IQ-TREE2 partitioned ML
JID3=$(sbatch --parsable --dependency=afterok:${JID2} "${SCRIPTS}/03_iqtree.sh")
echo "  Job 3 (IQ-TREE2):              ${JID3}"

echo ""
echo "Chain submitted. Monitor with:"
echo "  squeue -u nvillabo"
echo ""
echo "Final tree will be at:"
echo "  /dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v4/03_iqtree/snod_phylo_v4.treefile"
