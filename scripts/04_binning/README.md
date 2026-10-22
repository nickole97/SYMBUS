# Step 04 — Genome Binning & Refinement

Bins assembled contigs into MAGs using three independent binners (MetaBAT2, MaxBin2, CONCOCT), then integrates and refines them with DAS Tool to produce the best non-redundant set.

## Inputs
- `/work_meta/03_assembly/<sample>/final.contigs.fa`
- `/work_meta/04_mapping/<sample>/<sample>_depth.txt` (MetaBAT2, MaxBin2)
- `/work_meta/04_mapping/<sample>/<sample>.sorted.bam` (CONCOCT)

## Outputs
- `/work_meta/05_binning/<sample>/metabat/<sample>_metabat.*.fa`
- `/work_meta/05_binning/<sample>/maxbin/<sample>_maxbin.*.fasta`
- `/work_meta/05_binning/<sample>/concoct/bins/*.fa`
- `/work_meta/05_binning/<sample>/DASTool/<sample>_DASTool_DASTool_bins/*.fa` — final refined bins
- `/work_meta/05_binning/<sample>/DASTool/<sample>_DASTool_summary.txt`

## Software
- MetaBAT2 — tetranucleotide + coverage binning (`-m 1500`)
- MaxBin2 — tetranucleotide + abundance binning (`-min_contig_length 1500`)
- CONCOCT — composition + multi-sample coverage binning; contigs cut at 10 kb before clustering
- DAS Tool — bin refinement integrating all three binners (`--search_engine diamond`)
- CondaEnv: `binning_env` (all four tools)

## Usage
```bash
sbatch --array=1-113%5 4a_binning_metabat.sh
sbatch --array=1-113%5 4b_binning_maxbin.sh
sbatch --array=1-113%7 4c_binning_concoct.sh
sbatch --array=1-113   4d_refine_dastool.sh
```

## Notes
- Run scripts in order (4a → 4d); each step checks for existing output and skips completed samples
- DAS Tool requires at least one binner to have produced bins; samples with zero bins are skipped
- `Fasta_to_Contig2Bin.sh` converts bin FASTAs to TSV input for DAS Tool
