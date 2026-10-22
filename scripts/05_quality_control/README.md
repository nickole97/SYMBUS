# Step 05 — MAG Quality Assessment

Assesses completeness, contamination, and strain heterogeneity of DAS Tool bins using CheckM's lineage-specific marker gene workflow, then collects and merges results into a single table.

## Inputs
- `/work_meta/05_binning/<sample>/DASTool/<sample>_DASTool_DASTool_bins/*.fa`

## Outputs
- `/work_meta/05_binning/<sample>/CheckM/storage/bin_stats_ext.tsv` — per-bin detailed stats
- `/work_meta/05_binning/<sample>/CheckM/<sample>_checkm_qa.tsv` — formatted QA table (5b)
- `/work_meta/05_binning/checkm_qa_ALL.tsv` — combined table across all samples (5c)

## Software
- CheckM — lineage_wf (lineage-specific marker genes) + qa (format output)
- CondaEnv: `checkm_env` (with data at `$CHECKM_DATA_PATH`)

## Usage
```bash
# 5a: run per-sample (SLURM array)
sbatch --array=1-113 5a_checkm_ind.sh

# 5b: generate per-sample QA tables (local, no SLURM)
bash 5b_checkm_collect.sh

# 5c: merge all QA tables into one file (local, no SLURM)
bash 5c_checkm_combine.sh
```

## Notes
- `5a` runs `lineage_wf` on the DAS Tool bins folder; samples with no `.fa` files are skipped
- `5b` calls `checkm qa` with `-o 2 --tab_table` to produce a parseable TSV per sample
- `5c` adds a SAMPLE column to each table before concatenating into `checkm_qa_ALL.tsv`
- Quality thresholds for downstream selection (step 06): completeness ≥ 70%, contamination < 10%
