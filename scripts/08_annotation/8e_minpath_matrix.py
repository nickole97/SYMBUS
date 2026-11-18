#!/usr/bin/env python3
"""
8e_minpath_matrix.py
Build a MAG × pathway presence/absence matrix from MinPath .report files.

Output: minpath_matrix.tsv  (rows = MAGs, cols = KEGG pathways)
"""

import os
import re
import pandas as pd
from pathlib import Path

MINPATH_DIR = Path("/dfs10/hammert-lab/nvillabo/Ch2/work_meta/08_annotation/minpath")
OUT_FILE    = MINPATH_DIR / "minpath_matrix.tsv"

def parse_report(report_path):
    """Return dict {pathway_id: 1} for pathways called present by MinPath."""
    present = {}
    with open(report_path) as fh:
        for line in fh:
            # MinPath report format:
            # path 00010 kegg n/a  naive 1  minpath 1  fam0  105  fam-found  17  name  Glycolysis...
            m = re.match(r"path\s+(\S+)\s+\S+\s+\S+\s+naive\s+\d+\s+minpath\s+(\d)", line)
            if m and m.group(2) == "1":
                present[m.group(1)] = 1
    return present

records = {}
for mag_dir in sorted(MINPATH_DIR.iterdir()):
    report = mag_dir / "minpath.report"
    if not report.exists():
        print(f"[WARN] missing report: {mag_dir.name}")
        continue
    records[mag_dir.name] = parse_report(report)

if not records:
    raise SystemExit("No MinPath reports found — check MINPATH_DIR path.")

df = pd.DataFrame(records).T.fillna(0).astype(int)
df.index.name = "MAG"
df = df.sort_index()
df = df[sorted(df.columns)]

df.to_csv(OUT_FILE, sep="\t")
print(f"Matrix written: {df.shape[0]} MAGs × {df.shape[1]} pathways -> {OUT_FILE}")
