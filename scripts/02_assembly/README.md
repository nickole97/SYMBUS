# Step 02 — Metagenomic Assembly

Assembles non-host reads into contigs using MEGAHIT, discarding contigs shorter than 1000 bp.

## Inputs
- `/work_meta/preprocessing/02_host_removal/<sample>/*_R1_nonhost.fastq.gz`
- `/work_meta/preprocessing/02_host_removal/<sample>/*_R2_nonhost.fastq.gz`
- `/work_meta/preprocessing/02_host_removal/<sample>/*_singleton_nonhost.fastq.gz` (optional)

## Outputs
- `/work_meta/03_assembly/<sample>/final.contigs.fa` — assembled contigs ≥ 1000 bp

## Software
- MEGAHIT 1.2.9 — de novo metagenome assembler (`--min-contig-len 1000`)
- CondaEnv: loaded via `module load megahit/1.2.9`

## Usage
```bash
sbatch --array=0-112 2a_assembly_megahit.sh
```

## Notes
- Singleton reads (if present) are included via the `-r` flag
- Runtime is approximately 12–24 hours per sample depending on read depth
- Job skips samples where `final.contigs.fa` already exists
