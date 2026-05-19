#!/usr/bin/env python3
"""
Consolidates KofamScan ko_mapper.txt files into a MAG x KO count matrix.

Each cell = number of genes in that MAG assigned to that KO (0 if absent).
Mirrors the format of cazyme_matrix_all.tsv.

Usage:
    python 8c_kofamscan_matrix.py
    python 8c_kofamscan_matrix.py --kofamscan-dir /path/to/kofamscan --out /path/to/out.tsv
"""

import os
import sys
import argparse
from collections import defaultdict

KOFAMSCAN_DIR = "/dfs10/hammert-lab/nvillabo/Ch2/work_meta/08_annotation/kofamscan"
OUT_FILE      = "/dfs10/hammert-lab/nvillabo/Ch2/work_meta/08_annotation/ko_matrix_all.tsv"

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--kofamscan-dir", default=KOFAMSCAN_DIR)
    p.add_argument("--out", default=OUT_FILE)
    return p.parse_args()

def main():
    args = parse_args()
    kofamscan_dir = args.kofamscan_dir
    out_file      = args.out

    if not os.path.isdir(kofamscan_dir):
        sys.exit(f"ERROR: {kofamscan_dir} not found")

    mag_ko_counts = {}   # {MAG: {KO: count}}
    all_kos       = set()
    missing       = []

    for mag in sorted(os.listdir(kofamscan_dir)):
        mapper = os.path.join(kofamscan_dir, mag, "ko_mapper.txt")
        if not os.path.isfile(mapper):
            missing.append(mag)
            continue

        ko_counts = defaultdict(int)
        with open(mapper) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    ko = parts[1]
                    ko_counts[ko] += 1
                    all_kos.add(ko)

        mag_ko_counts[mag] = dict(ko_counts)

    if missing:
        print(f"WARNING: {len(missing)} MAG(s) missing ko_mapper.txt — skipped:", file=sys.stderr)
        for m in missing:
            print(f"  {m}", file=sys.stderr)

    mags = sorted(mag_ko_counts.keys())
    kos  = sorted(all_kos)

    print(f"Building matrix: {len(mags)} MAGs x {len(kos)} KOs -> {out_file}", file=sys.stderr)

    with open(out_file, "w") as fh:
        fh.write("MAG\t" + "\t".join(kos) + "\n")
        for mag in mags:
            counts = mag_ko_counts[mag]
            row = [mag] + [str(counts.get(ko, 0)) for ko in kos]
            fh.write("\t".join(row) + "\n")

    print(f"Done: {len(mags)} MAGs, {len(kos)} unique KOs", file=sys.stderr)

if __name__ == "__main__":
    main()
