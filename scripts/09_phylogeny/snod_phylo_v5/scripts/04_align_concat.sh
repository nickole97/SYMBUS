#!/usr/bin/env bash
#SBATCH --job-name=align_snodv5
#SBATCH --output=logs/04_align_%j.out
#SBATCH --error=logs/04_align_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab

source ~/.bashrc
conda activate phylo
set -euo pipefail

ROOT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v5"
OG_DIR="${ROOT}/04_core_OGs"
ALN_DIR="${ROOT}/05_alignments"
mkdir -p "${ALN_DIR}/raw" "${ALN_DIR}/trimmed"

N_CORE=$(wc -l < "${OG_DIR}/core_OG_ids.txt")
echo "Aligning ${N_CORE} core orthogroups..."

while IFS= read -r OG; do
    # Skip OGs already aligned and trimmed (rerun after the concatenation fix).
    if [[ -s "${ALN_DIR}/trimmed/${OG}.trimmed.aln" ]]; then continue; fi
    mafft \
        --auto \
        --thread ${SLURM_CPUS_PER_TASK} \
        --quiet \
        "${OG_DIR}/${OG}.faa" \
        > "${ALN_DIR}/raw/${OG}.aln"

    trimal \
        -in  "${ALN_DIR}/raw/${OG}.aln" \
        -out "${ALN_DIR}/trimmed/${OG}.trimmed.aln" \
        -automated1
done < "${OG_DIR}/core_OG_ids.txt"

echo "Alignments done. Concatenating..."

python3 - <<'PYEOF'
import os
from pathlib import Path
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict

root    = Path("/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/snod_phylo_v5")
aln_dir = root / "05_alignments/trimmed"
out_dir = root / "05_alignments"

all_alns = sorted(aln_dir.glob("*.trimmed.aln"))
print(f"Found {len(all_alns)} trimmed alignments")

# With --min_taxa < all, an OG can be absent from some genomes. Every taxon
# must still get a block for every partition, so absent taxa are padded with
# gaps; otherwise sequence lengths differ and IQ-TREE refuses the matrix
# ("contains not enough characters", first v5 run 2026-11-25).
taxa = set()
for aln_file in all_alns:
    for rec in AlignIO.read(aln_file, "fasta"):
        taxa.add(rec.id)
taxa = sorted(taxa)
print(f"Taxa across alignments: {len(taxa)}")

supermatrix     = {t: "" for t in taxa}
partition_lines = []
pos = 1
n_padded = 0

for aln_file in all_alns:
    aln     = AlignIO.read(aln_file, "fasta")
    og      = aln_file.stem.replace(".trimmed", "")
    aln_len = aln.get_alignment_length()
    partition_lines.append(f"LG, {og} = {pos}-{pos + aln_len - 1}")
    pos += aln_len
    present = {rec.id: str(rec.seq) for rec in aln}
    for t in taxa:
        if t in present:
            supermatrix[t] += present[t]
        else:
            supermatrix[t] += "-" * aln_len
            n_padded += 1
print(f"Padded (taxon, OG) cells: {n_padded}")
lengths = {len(v) for v in supermatrix.values()}
assert len(lengths) == 1, f"unequal lengths after padding: {lengths}"

concat_fasta = out_dir / "supermatrix.faa"
records = [SeqRecord(Seq(seq), id=name, description="")
           for name, seq in supermatrix.items()]
SeqIO.write(records, concat_fasta, "fasta")
print(f"Supermatrix: {concat_fasta}")
print(f"  Taxa: {len(records)}, Length: {len(records[0].seq)} aa")

partition_file = out_dir / "partitions.txt"
with open(partition_file, "w") as fh:
    fh.write("\n".join(partition_lines) + "\n")
print(f"Partition file: {partition_file}")
PYEOF

echo "Concatenation complete."
ls -lh "${ROOT}/05_alignments/supermatrix.faa" "${ROOT}/05_alignments/partitions.txt"
