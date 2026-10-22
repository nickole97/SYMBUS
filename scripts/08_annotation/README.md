# Step 08 — Functional & Taxonomic Annotation

Assigns taxonomy to MAGs with GTDB-Tk (8a) and annotates genes with KEGG Orthology using KofamScan (8b), then builds a MAG × KO count matrix (8c).

## Inputs
- `/work_meta/06_MAGs_selected/*.fa` — selected MAGs (8a, 8b setup)
- `/work_meta/08_annotation/dbcan_all/<MAG>/uniInput.faa` — protein predictions from dbCAN/Prodigal (8b)
- `/dfs10/hammert-lab/nvillabo/databases/kofamscan/` — KEGG HMM profiles (downloaded by 8b_setup)

## Outputs
- `/work_meta/08_gtdbtk/gtdbtk.bac120.summary.tsv` — bacterial taxonomy
- `/work_meta/08_gtdbtk/gtdbtk.ar53.summary.tsv` — archaeal taxonomy
- `/work_meta/08_annotation/kofamscan/<MAG>/ko_mapper.txt` — gene→KO assignments per MAG
- `/work_meta/08_annotation/ko_matrix_all.tsv` — MAG × KO count matrix (8c)

## Software
- GTDB-Tk 2.5.2 — `classify_wf`; phylogenetic placement against GTDB r226
- KofamScan (`exec_annotation`) — HMM-based KO annotation (`-f mapper`)
- Python 3 (stdlib + pandas) — 8c matrix consolidation script
- CondaEnv: `gtdbtk-2.5.2` (8a); `kofamscan_env` (8b, installed by 8b_setup)

## Usage
```bash
# Taxonomy
sbatch 8a_taxonomy_gtdbtk.sh

# KofamScan: one-time setup, then array
sbatch 8b_kofamscan_setup.sh
sbatch --array=1-121 8b_kofamscan.sh
# For timed-out MAGs (11 MAGs, 24 h limit):
sbatch --array=1-11  8b_kofamscan_rerun.sh

# Build KO matrix
python 8c_kofamscan_matrix.py
```

## Notes
- `8a` requires ~256 GB RAM for skani ANI screening against all 143k GTDB r226 reference genomes
- `8b_rerun` uses `kofamscan_missing_list.txt` instead of the main list; it was created for 11 MAGs that exceeded the original 12 h time limit
- `8c` accepts `--kofamscan-dir` and `--out` flags to override default paths
- KofamScan uses Ruby; the script retries up to 3× on NFS cache-miss startup failures
