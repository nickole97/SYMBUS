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

with open(args.orthogroups) as fh:
    header = fh.readline().strip().split("\t")
    taxa = header[1:]
    n_taxa = len(taxa)

min_taxa = args.min_taxa if args.min_taxa is not None else n_taxa
print(f"Taxa in OrthoFinder table: {n_taxa}")
print(f"Requiring presence in:      {min_taxa} (all={n_taxa})")

core_OGs = []
multi_copy = 0
missing = 0

with open(args.orthogroups) as fh:
    fh.readline()  # skip header
    for line in fh:
        parts = line.strip().split("\t")
        og_id = parts[0]
        genes = parts[1:]

        present = [g for g in genes if g.strip()]
        if len(present) < min_taxa:
            missing += 1
            continue

        if any("," in g for g in present if g.strip()):
            multi_copy += 1
            continue

        core_OGs.append(og_id)

total = len(core_OGs) + multi_copy + missing
print(f"\nTotal OGs parsed:       {total}")
print(f"  Multi-copy (skip):    {multi_copy}")
print(f"  Missing taxa (skip):  {missing}")
print(f"  Single-copy core:     {len(core_OGs)}")

faa_dir = Path(args.faa_dir)
not_found = []
written = 0
for og in core_OGs:
    src = faa_dir / f"{og}.fa"
    if not src.exists():
        not_found.append(og)
        continue
    records = []
    for rec in SeqIO.parse(src, "fasta"):
        # OrthoFinder gene IDs: genome_id_geneN
        # Strip the trailing _N (gene number) to recover the genome_id.
        # This handles genome IDs that themselves contain underscores.
        parts = rec.id.rsplit("_", 1)
        genome_id = parts[0] if len(parts) == 2 and parts[1].isdigit() else rec.id
        rec.id = genome_id
        rec.description = ""
        records.append(rec)
    SeqIO.write(records, outdir / f"{og}.faa", "fasta")
    written += 1

print(f"\nWrote {written} OG .faa files to: {outdir}")
if not_found:
    print(f"WARNING — OG .fa files not found: {not_found}")

og_list_file = outdir / "core_OG_ids.txt"
with open(og_list_file, "w") as fh:
    fh.write("\n".join(core_OGs) + "\n")
print(f"OG ID list: {og_list_file}")
