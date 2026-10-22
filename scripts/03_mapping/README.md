# Step 03 — Read Mapping to Contigs

Maps non-host reads back to their sample's assembled contigs to generate coverage depth information required for binning.

## Inputs
- `/work_meta/preprocessing/02_host_removal/<sample>/*_R1/R2_nonhost.fastq.gz`
- `/work_meta/03_assembly/<sample>/final.contigs.fa`

## Outputs
- `/work_meta/04_mapping/<sample>/<sample>.sorted.bam` — sorted BAM (+ `.bai` index)
- `/work_meta/04_mapping/<sample>/<sample>_depth.txt` — per-contig coverage depth table

## Software
- Minimap2 2.28 — short-read mapping (`-ax sr`) against per-sample contig index
- Samtools 1.15.1 — SAM→BAM conversion, sorting, and indexing
- jgi_summarize_bam_contig_depths (MetaBAT2) — generates the depth table used by binners
- CondaEnv: `quast_env` (provides the jgi_summarize tool)

## Usage
```bash
sbatch --array=1-113%5 3a_mapping_minimap2.sh
```

## Notes
- The Minimap2 index (`.mmi`) is built per sample and reused if it already exists
- The `%5` throttle limits concurrent jobs to 5 to avoid I/O contention
- The depth file format is the same as MetaBAT2 expects for binning input
