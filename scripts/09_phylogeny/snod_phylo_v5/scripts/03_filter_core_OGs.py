#!/usr/bin/env python3
"""
03_filter_core_OGs.py
Filter OrthoFinder output to single-copy core orthogroups
present in ALL genomes (or --min_taxa).

Usage:
    python3 03_filter_core_OGs.py \
        --orthogroups <OrthoFinder results>/Orthogroups/Orthogroups.tsv \
        --faa_dir     <OrthoFinder results>/Orthogroup_Sequences/ \
        --outdir      ~/snod_phylo_v2/05_core_OGs \
        --min_taxa    N    # N = total number of genomes (default: all)
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--orthogroups", required=True)
parser.add_argument("--faa_dir", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--min_taxa", type=int, default=None,
                    help="Min taxa required; default = all taxa in the table")
args = parser.parse_args()

outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

# Column headers in Orthogroups.tsv ARE the genome IDs (input .faa basenames)
with open(args.orthogroups) as fh:
    header = fh.readline().strip().split("\t")
taxa = header[1:]
n_taxa = len(taxa)

min_taxa = args.min_taxa if args.min_taxa is not None else n_taxa
print(f"Taxa in OrthoFinder table: {n_taxa}")
print(f"Requiring presence in:      {min_taxa} (all={n_taxa})")

core_OGs = []   # list of (og_id, {genome_id: gene_id})
multi_copy = 0
missing = 0

with open(args.orthogroups) as fh:
    fh.readline()  # skip header
    for line in fh:
        parts = line.strip().split("\t")
        og_id = parts[0]
        genes = parts[1:]
        # TSV may omit trailing tabs for absent taxa
        while len(genes) < n_taxa:
            genes.append("")

        present = [i for i, g in enumerate(genes) if g.strip()]
        if len(present) < min_taxa:
            missing += 1
            continue

        if any("," in genes[i] for i in present):
            multi_copy += 1
            continue

        genome_to_gene = {taxa[i]: genes[i].strip() for i in present}
        core_OGs.append((og_id, genome_to_gene))

total = len(core_OGs) + multi_copy + missing
print(f"\nTotal OGs parsed:       {total}")
print(f"  Multi-copy (skip):    {multi_copy}")
print(f"  Missing taxa (skip):  {missing}")
print(f"  Single-copy core:     {len(core_OGs)}")

faa_dir = Path(args.faa_dir)
not_found = []
written = 0

for og_id, genome_to_gene in core_OGs:
    src = faa_dir / f"{og_id}.fa"
    if not src.exists():
        not_found.append(og_id)
        continue

    # Index all sequences by their record ID (full gene header after sed)
    seq_index = SeqIO.to_dict(SeqIO.parse(src, "fasta"))

    records = []
    for genome_id, gene_id in genome_to_gene.items():
        if gene_id not in seq_index:
            print(f"WARNING: '{gene_id}' not in {src.name}", file=sys.stderr)
            continue
        rec = seq_index[gene_id]
        rec.id = genome_id
        rec.description = ""
        records.append(rec)

    SeqIO.write(records, outdir / f"{og_id}.faa", "fasta")
    written += 1

print(f"\nWrote {written} OG .faa files to: {outdir}")
if not_found:
    print(f"WARNING — OG .fa files not found: {not_found}")

og_list_file = outdir / "core_OG_ids.txt"
with open(og_list_file, "w") as fh:
    fh.write("\n".join(og for og, _ in core_OGs) + "\n")
print(f"OG ID list: {og_list_file}")
