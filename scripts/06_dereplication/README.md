# Step 06 — MAG Selection & Dereplication

Copies quality-passing MAGs to a central directory, then dereplicates them at 99% ANI (species level) with dRep to produce one representative genome per species.

## Inputs
- `/work_meta/MagsToUse.tsv` — two-column TSV (SAMPLE, BINID) listing MAGs that pass quality thresholds
- `/work_meta/05_binning/<sample>/DASTool/<sample>_DASTool_DASTool_bins/<binid>.fa`

## Outputs
- `/work_meta/06_MAGs_selected/<sample>__<binid>.fa` — selected MAGs, one file each (6a)
- `/work_meta/07_dRep_99_v3/dereplicated_genomes/*.fa` — one representative per species (6b)
- `/work_meta/07_dRep_99_v3/data_tables/` — clustering and quality tables
- `/work_meta/07_dRep_99_v3/figures/` — dendrograms and summary plots

## Software
- dRep — dereplication (`-pa 0.9` primary MASH clustering, `-sa 0.99` secondary fastANI, `-clusterAlg single`)
- CondaEnv: `drep`

## Usage
```bash
# 6a: run locally — copies FASTAs to the selection directory
bash 6a_select_mags.sh

# 6b: submit dereplication job
sbatch 6b_drep.sh
```

## Notes
- `6a` searches multiple candidate paths for each bin (metabat, maxbin, concoct, DASTool) in priority order; missing bins are flagged with a warning
- dRep thresholds: primary 90% (MASH pre-filter), secondary 99% ANI (fastANI); single-linkage clustering
- Estimated runtime for 6b: 6–12 hours depending on number of input MAGs
