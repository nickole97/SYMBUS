"""
Consolidate dbCAN5 overview.tsv files from multiple MAGs into a single CAZyme count matrix.
Filters for hits with #ofTools >= 2 (high confidence).

Usage:
    python consolidate_cazymes.py \
        --input_dir /path/to/dbcan/output/ \
        --output /path/to/cazyme_matrix.tsv

Output:
    TSV with MAG names as rows, CAZyme families as columns, gene counts as values.
"""

import os
import argparse
import pandas as pd
from collections import defaultdict

def parse_overview(filepath, min_tools=2):
    """
    Parse a dbCAN overview.tsv file and return a dict of {CAZyme_family: count}.
    Only counts genes with #ofTools >= min_tools.
    """
    counts = defaultdict(int)
    
    try:
        df = pd.read_csv(filepath, sep='\t')
    except Exception as e:
        print(f"  WARNING: Could not read {filepath}: {e}")
        return counts
    
    # Column names in dbCAN5 overview.tsv
    # Gene ID | EC# | HMMER | dbCAN_sub | DIAMOND | #ofTools | Substrate
    tools_col = '#ofTools'
    hmmer_col = 'dbCAN_hmm'
    
    if tools_col not in df.columns:
        print(f"  WARNING: '{tools_col}' column not found in {filepath}")
        print(f"  Available columns: {list(df.columns)}")
        return counts
    
    # Filter by number of tools
    df_filtered = df[df[tools_col] >= min_tools].copy()
    
    if hmmer_col not in df_filtered.columns:
        print(f"  WARNING: '{hmmer_col}' column not found in {filepath}")
        return counts
    
    # Extract CAZyme families from HMMER column
    # HMMER entries look like: GH43_4(2-398)+GH43_27(400-600) or just GH43
    for entry in df_filtered[hmmer_col].dropna():
        if entry == '-':
            continue
        # Split multiple domains (separated by '+')
        domains = entry.split('+')
        for domain in domains:
            # Remove position info in parentheses: GH43_4(2-398) -> GH43_4
            family = domain.split('(')[0].strip()
            if family and family != '-':
                # Optionally collapse subfamilies: GH43_4 -> GH43
                # We keep full subfamily for now (more informative)
                counts[family] += 1
    
    return counts


def main():
    parser = argparse.ArgumentParser(description='Consolidate dbCAN CAZyme annotations across MAGs')
    parser.add_argument('--input_dir', required=True,
                        help='Directory containing per-MAG subdirectories with overview.tsv files')
    parser.add_argument('--output', required=True,
                        help='Output TSV file path for the consolidated matrix')
    parser.add_argument('--min_tools', type=int, default=2,
                        help='Minimum number of tools that must agree (default: 2)')
    parser.add_argument('--collapse_subfamilies', action='store_true',
                        help='Collapse subfamilies (e.g. GH43_4 -> GH43)')
    args = parser.parse_args()

    input_dir = args.input_dir.rstrip('/')
    print(f"\nScanning: {input_dir}")
    print(f"Min tools filter: >= {args.min_tools}")
    print(f"Collapse subfamilies: {args.collapse_subfamilies}\n")

    # Find all overview.tsv files
    mag_counts = {}
    
    subdirs = sorted([d for d in os.listdir(input_dir) 
                      if os.path.isdir(os.path.join(input_dir, d))])
    
    print(f"Found {len(subdirs)} MAG directories\n")
    
    for subdir in subdirs:
        overview_path = os.path.join(input_dir, subdir, 'overview.tsv')
        
        if not os.path.exists(overview_path):
            print(f"  SKIP: No overview.tsv in {subdir}")
            continue
        
        counts = parse_overview(overview_path, min_tools=args.min_tools)
        
        if args.collapse_subfamilies:
            collapsed = defaultdict(int)
            for family, count in counts.items():
                base = family.split('_')[0]
                collapsed[base] += count
            counts = collapsed
        
        mag_counts[subdir] = counts
        total = sum(counts.values())
        print(f"  {subdir}: {total} CAZyme genes ({len(counts)} unique families)")
    
    print(f"\nBuilding matrix for {len(mag_counts)} MAGs...")
    
    # Build dataframe: rows = MAGs, columns = CAZyme families
    df = pd.DataFrame(mag_counts).T.fillna(0).astype(int)
    df.index.name = 'MAG'
    
    # Sort columns by CAZyme class then number
    def sort_key(col):
        import re
        m = re.match(r'([A-Z]+)(\d+)', col)
        if m:
            return (m.group(1), int(m.group(2)))
        return (col, 0)
    
    df = df[sorted(df.columns, key=sort_key)]
    
    # Save
    df.to_csv(args.output, sep='\t')
    print(f"\nSaved matrix: {args.output}")
    print(f"  Shape: {df.shape[0]} MAGs x {df.shape[1]} CAZyme families")
    print(f"\nTop 10 most abundant families:")
    print(df.sum().sort_values(ascending=False).head(10).to_string())


if __name__ == '__main__':
    main()
