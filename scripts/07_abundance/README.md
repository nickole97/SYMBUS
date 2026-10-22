# Step 07 — MAG Abundance Quantification

Maps non-host reads from each sample to the full set of selected MAGs to calculate per-genome abundance metrics (mean coverage, covered fraction, TPM).

## Inputs
- `/work_meta/preprocessing/02_host_removal/<sample>/*_R1/R2_nonhost.fastq.gz`
- `/work_meta/06_MAGs_selected/*.fa` — all selected MAGs as reference

## Outputs
- `/work_meta/07_abundance/<sample>_coverm_genome.tsv` — per-MAG abundance table for each sample

## Software
- CoverM (genome mode) — read mapping and coverage calculation
  - Mapper: `minimap2-sr`
  - Metrics: `mean`, `covered_fraction`, `tpm`
  - Filters: `--min-read-percent-identity 95`, `--min-read-aligned-percent 75`
- CondaEnv: `coverm_env`

## Usage
```bash
sbatch --array=1-113 7a_coverm_genome.sh
```

## Notes
- This step quantifies how much each sample's reads map to each dereplicated MAG, producing a cross-sample abundance matrix when the per-sample TSVs are merged downstream
- The script notes that this approach may be replaced by inStrain (step 10) for strain-level profiling
- Jobs skip samples where the output TSV already exists
