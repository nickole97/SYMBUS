#!/usr/bin/env python3
"""
Concatenate trimmed OG alignments into a supermatrix.
Taxa missing from an OG are filled with gap characters (-).
Writes supermatrix.faa and partitions.txt.
"""
import argparse
from pathlib import Path
from collections import defaultdict, OrderedDict
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--trimmed_dir", required=True)
    p.add_argument("--outdir",      required=True)
    return p.parse_args()

def main():
    args = parse_args()
    trimmed_dir = Path(args.trimmed_dir)
    outdir      = Path(args.outdir)

    aln_files = sorted(trimmed_dir.glob("*.trimmed.aln"))
    print(f"Found {len(aln_files)} trimmed alignments")

    # First pass: collect all taxa names across all OGs
    all_taxa = set()
    og_data  = OrderedDict()   # og_id → {taxon → sequence}
    og_lens  = {}

    for aln_file in aln_files:
        og = aln_file.stem.replace(".trimmed", "")
        aln = AlignIO.read(aln_file, "fasta")
        aln_len = aln.get_alignment_length()
        og_lens[og] = aln_len
        og_data[og] = {}
        for rec in aln:
            # Strip OrthoFinder-style suffixes so taxon names are consistent
            # Ingroup: header looks like GenomeID_1234; Outgroup: just the name
            taxon = rec.id
            og_data[og][taxon] = str(rec.seq)
            all_taxa.add(taxon)

    all_taxa = sorted(all_taxa)
    print(f"Total taxa (ingroup + outgroup): {len(all_taxa)}")

    # Second pass: build supermatrix with gap-filling
    supermatrix = {t: "" for t in all_taxa}
    partition_lines = []
    pos = 1

    for og, seqs in og_data.items():
        aln_len = og_lens[og]
        partition_lines.append(f"LG, {og} = {pos}-{pos + aln_len - 1}")
        pos += aln_len
        for taxon in all_taxa:
            if taxon in seqs:
                supermatrix[taxon] += seqs[taxon]
            else:
                supermatrix[taxon] += "-" * aln_len

    # Write supermatrix FASTA
    out_fasta = outdir / "supermatrix.faa"
    records = [SeqRecord(Seq(supermatrix[t]), id=t, description="")
               for t in all_taxa]
    SeqIO.write(records, out_fasta, "fasta")
    print(f"Supermatrix written: {out_fasta}")
    print(f"  Taxa:   {len(records)}")
    if records:
        print(f"  Length: {len(records[0].seq)} aa")

    # Write partition file
    out_part = outdir / "partitions.txt"
    with open(out_part, "w") as fh:
        fh.write("\n".join(partition_lines) + "\n")
    print(f"Partition file written: {out_part}")
    print(f"  Partitions: {len(partition_lines)}")

if __name__ == "__main__":
    main()
