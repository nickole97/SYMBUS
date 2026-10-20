#!/usr/bin/env bash
#SBATCH --job-name=align_snodv2
#SBATCH --output=logs/04_align_%j.out
#SBATCH --error=logs/04_align_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=standard

# ─────────────────────────────────────────────────────────────────────────────
# Step 4: Align each core orthogroup (MAFFT), trim (trimAl),
#         then concatenate into a protein supermatrix
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail
source ~/.bashrc
conda activate mafft-env 2>/dev/null || conda activate base

ROOT="$HOME/snod_phylo_v2"
OG_DIR="${ROOT}/05_core_OGs"
ALN_DIR="${ROOT}/06_alignments"
mkdir -p "${ALN_DIR}/raw" "${ALN_DIR}/trimmed"

N_CORE=$(wc -l < "${OG_DIR}/core_OG_ids.txt")
echo "Aligning ${N_CORE} core orthogroups..."

while IFS= read -r OG; do
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

root = Path(os.environ["HOME"]) / "snod_phylo_v2"
aln_dir = root / "06_alignments/trimmed"
out_dir  = root / "06_alignments"

all_alns = sorted(aln_dir.glob("*.trimmed.aln"))
print(f"Found {len(all_alns)} trimmed alignments")

supermatrix = defaultdict(str)
partition_lines = []
pos = 1

for aln_file in all_alns:
    aln = AlignIO.read(aln_file, "fasta")
    og  = aln_file.stem.replace(".trimmed", "")
    aln_len = aln.get_alignment_length()
    partition_lines.append(f"LG, {og} = {pos}-{pos + aln_len - 1}")
    pos += aln_len
    for rec in aln:
        supermatrix[rec.id] += str(rec.seq)

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
ls -lh "${ROOT}/06_alignments/supermatrix.faa" "${ROOT}/06_alignments/partitions.txt"
