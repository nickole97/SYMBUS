# Step 10 — Within-Sample Strain Profiling (inStrain)

Profiles nucleotide-level diversity and calculates population ANI (popANI) across all samples for each dereplicated MAG using inStrain, enabling strain-level comparisons between bee hosts.

## Inputs
- `/work_meta/07_dRep_95_correct/dereplicated_genomes/*.fa` — dereplicated MAGs
- `/work_meta/preprocessing/clean_reads/*_R1/R2_nonhost.fastq.gz` — non-host reads per sample

## Outputs
- `10_inStrain/dereplicated_genomes.fasta` — concatenated MAGs with renamed scaffolds (10b)
- `10_inStrain/dereplicated_genomes.stb` — scaffold-to-bin file (10b)
- `10_inStrain/dereplicated_genomes.genes.fna/.faa` — Prodigal gene predictions (10a)
- `10_inStrain/10b_mapping/<sample>.sorted.bam` — per-sample BAMs mapped to concat reference (10c)
- `10_inStrain/10c_profiles/<sample>.IS/` — per-sample inStrain profiles (10d)
- `10_inStrain/10d_compare/all_samples/` — pairwise popANI comparisons (10e)

## Software
- Prodigal — gene prediction (meta mode) over concatenated reference (10a)
- Minimap2 — short-read mapping (`-ax sr`) to concatenated reference (10c)
- Samtools — BAM sorting and indexing (10c)
- inStrain — `profile` (per-sample; `-l 0.95 -c 1 -f 0.05`) and `compare` (all-vs-all, `--database_mode`) (10d, 10e)
- CondaEnv: `instrain_env` (10c–10e); `binning_env` (10a, 10b)

## Usage
```bash
# Run steps in order
sbatch 10b_prepare_refs.sh   # builds concat FASTA and STB
sbatch 10a_prodigal.sh       # gene calling (depends on 10b output)
sbatch --array=1-113 10c_mapping.sh
sbatch --array=1-113 10d_instrain_profile.sh
sbatch 10e_instrain_compare.sh
```

## Notes
- Scaffold headers are renamed `BINNAME__original_header` in the concat FASTA to avoid duplicate names; the STB file maps these back to bin names
- inStrain profile thresholds: breadth filter `-l 0.95`, min coverage `-c 1`, min allele frequency `-f 0.05`
- `10e` requires all profiles to exist; run after all `10d` array jobs complete
- Scripts are written in Spanish (comments/variable names) — funcionality is identical
