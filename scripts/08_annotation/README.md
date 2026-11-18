# Step 08 — Functional & Taxonomic Annotation

Assigns taxonomy to MAGs with GTDB-Tk (8a), annotates genes with KEGG Orthology using KofamScan (8b) and builds a MAG × KO count matrix (8c), then reconstructs KEGG pathways per MAG with MinPath (8d) and builds a MAG × pathway matrix (8e).

## Inputs
- `/work_meta/06_MAGs_selected/*.fa` — selected MAGs (8a, 8b setup)
- `/work_meta/08_annotation/dbcan_all/<MAG>/uniInput.faa` — protein predictions from dbCAN/Prodigal (8b)
- `/dfs10/hammert-lab/nvillabo/databases/kofamscan/` — KEGG HMM profiles (downloaded by 8b_setup)
- `/work_meta/08_annotation/kofamscan/<MAG>/ko_mapper.txt` — KO assignments, consumed by 8d
- `/work_meta/08_annotation/minpath_mag_list.txt` — 121 MAGs, one per line (8d array index)

## Outputs
- `/work_meta/08_gtdbtk/gtdbtk.bac120.summary.tsv` — bacterial taxonomy
- `/work_meta/08_gtdbtk/gtdbtk.ar53.summary.tsv` — archaeal taxonomy
- `/work_meta/08_annotation/kofamscan/<MAG>/ko_mapper.txt` — gene→KO assignments per MAG
- `/work_meta/08_annotation/ko_matrix_all.tsv` — MAG × KO count matrix (8c)
- `/work_meta/08_annotation/minpath/<MAG>/minpath.report` — pathways per MAG; `minpath 1` = present after parsimony, `naive 1` = present by any gene hit (8d)
- `/work_meta/08_annotation/minpath/<MAG>/minpath.details` — KOs supporting each pathway (8d)
- `/work_meta/08_annotation/minpath/minpath_matrix.tsv` — MAG × pathway presence/absence, 121 MAGs × 140 pathways (8e)

## Software
- GTDB-Tk 2.5.2 — `classify_wf`; phylogenetic placement against GTDB r226
- KofamScan (`exec_annotation`) — HMM-based KO annotation (`-f mapper`)
- MinPath — parsimony-based pathway reconstruction (integer linear programming); needs the GLPK solver
- Python 3 (stdlib + pandas) — 8c and 8e matrix consolidation scripts
- CondaEnv: `gtdbtk-2.5.2` (8a); `kofamscan_env` (8b, installed by 8b_setup); 8d uses the system `python3`, no env

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

# MinPath pathway reconstruction, then pathway matrix
sbatch --array=1-121 8d_minpath.sh
python 8e_minpath_matrix.py
```

## MinPath installation (manual — no setup script)

`8d` expects MinPath at `/dfs10/hammert-lab/nvillabo/tools/MinPath/MinPath.py`. It was installed by hand on 2026-06-07 and **no script in this repo reproduces it**. The bioconda `minpath` package was tried first and did not produce a working `MinPath.py`, so it was abandoned. To reinstall, from an interactive session (`inter`, never the login node):

```bash
mkdir -p /dfs10/hammert-lab/nvillabo/tools && cd /dfs10/hammert-lab/nvillabo/tools
git clone https://github.com/mgtools/MinPath    # pinned at 46d3e81
cd MinPath/glpk-4.6 && ./configure && make      # builds the glpsol solver MinPath calls
```

## Notes
- `8a` requires ~256 GB RAM for skani ANI screening against all 143k GTDB r226 reference genomes
- `8b_rerun` uses `kofamscan_missing_list.txt` instead of the main list; it was created for 11 MAGs that exceeded the original 12 h time limit
- `8c` accepts `--kofamscan-dir` and `--out` flags to override default paths
- KofamScan uses Ruby; the script retries up to 3× on NFS cache-miss startup failures
- `8d` must `cd` into the per-MAG output directory before running MinPath: MinPath writes temp files (`test.mps`, `test.mps.LPout`) to the current working directory, so parallel array tasks sharing a cwd would overwrite each other
- `8d` skips MAGs with zero assigned KOs rather than failing
- MinPath's parsimony is deliberately conservative: it can drop genuinely present pathways whose enzymes are fully shared with other pathways, and coverage is weaker for novel lineages poorly represented in KEGG
