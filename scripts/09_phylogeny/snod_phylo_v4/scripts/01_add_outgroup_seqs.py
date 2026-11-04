#!/usr/bin/env python3
"""
Parse BLAST results (from 00_blast_outgroups.sh) and add best-hit outgroup
sequences to each core OG FAA. Output FAAs contain ingroup sequences + any
outgroup sequence with a valid hit (evalue < threshold).
Outgroups without a hit for a given OG are omitted; the concat step fills gaps.
"""
import argparse
import os
from pathlib import Path
from Bio import SeqIO

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--v2_og_dir",    required=True, help="v2 core OG FAA directory (411 OGs)")
    p.add_argument("--blast_dir",    required=True, help="00_blast directory with {outgroup}_hits/")
    p.add_argument("--out_faa_dir",  required=True, help="v3 02_prodigal dir with outgroup FAAs")
    p.add_argument("--outdir",       required=True, help="Output directory for augmented OG FAAs")
    p.add_argument("--evalue",       type=float, default=1e-5)
    return p.parse_args()

def load_faa_index(faa_path):
    return SeqIO.to_dict(SeqIO.parse(faa_path, "fasta"))

def best_hit(tsv_path, evalue_thresh):
    """Return (sseqid, evalue) of the best BLAST hit, or None if no hit passes threshold."""
    best = None
    best_score = -1
    with open(tsv_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            sseqid  = parts[1]
            evalue  = float(parts[4])
            bitscore = float(parts[5])
            if evalue <= evalue_thresh and bitscore > best_score:
                best_score = bitscore
                best = (sseqid, evalue)
    return best

def main():
    args = parse_args()
    v2_og_dir   = Path(args.v2_og_dir)
    blast_dir   = Path(args.blast_dir)
    out_faa_dir = Path(args.out_faa_dir)
    outdir      = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    og_list = sorted(v2_og_dir.glob("*.faa"))
    print(f"Input OGs: {len(og_list)}")

    # Discover outgroups from blast_dir subdirectory names
    outgroups = [d.name.replace("_hits", "") for d in blast_dir.iterdir()
                 if d.is_dir() and d.name.endswith("_hits")]
    print(f"Outgroups: {outgroups}")

    # Load outgroup FAA indexes (lazy)
    og_indexes = {}
    for og_name in outgroups:
        faa = out_faa_dir / f"{og_name}.faa"
        if faa.exists():
            og_indexes[og_name] = load_faa_index(faa)
            print(f"  Loaded {og_name}: {len(og_indexes[og_name])} proteins")
        else:
            print(f"  WARNING: {faa} not found — skipping outgroup {og_name}")

    hit_counts = {og: 0 for og in outgroups}
    total_ogs = len(og_list)

    for faa_path in og_list:
        og_id = faa_path.stem
        # Copy ingroup sequences
        ingroup_records = list(SeqIO.parse(faa_path, "fasta"))
        augmented = list(ingroup_records)

        for og_name in outgroups:
            if og_name not in og_indexes:
                continue
            tsv = blast_dir / f"{og_name}_hits" / f"{og_id}.tsv"
            if not tsv.exists() or tsv.stat().st_size == 0:
                continue
            hit = best_hit(tsv, args.evalue)
            if hit is None:
                continue
            sseqid, evalue = hit
            if sseqid not in og_indexes[og_name]:
                print(f"  WARNING: {sseqid} not found in {og_name} FAA")
                continue
            rec = og_indexes[og_name][sseqid]
            # Rename to outgroup name so concat can identify taxon
            new_rec = rec.__class__(rec.seq, id=og_name, description="")
            augmented.append(new_rec)
            hit_counts[og_name] += 1

        out_file = outdir / f"{og_id}.faa"
        SeqIO.write(augmented, out_file, "fasta")

    print(f"\nAugmented OGs written: {total_ogs}")
    for og_name, n in hit_counts.items():
        print(f"  {og_name}: {n}/{total_ogs} OGs with hit ({100*n/total_ogs:.1f}%)")

if __name__ == "__main__":
    main()
