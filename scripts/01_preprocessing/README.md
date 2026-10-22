# Step 01 — Preprocessing & Host Removal

Quality-trim raw paired-end reads with BBDuk/BBMerge, then remove host (*Bombus impatiens*) reads by mapping against the reference genome and retaining only unmapped (microbial) reads.

## Inputs
- `/Ch2/raw_01/raw_paired/<sample>/<sample>_{1,2}.fq.gz` — raw paired-end FASTQ files

## Outputs
**1a (QC):** `/work_meta/preprocessing/01_preprocessing/<sample>/`
- `*_R1_clean.fastq.gz`, `*_R2_clean.fastq.gz` — adapter-trimmed, quality-filtered reads
- `*_merged.fastq.gz`, `*_R1/R2_unmerged.fastq.gz` — BBMerge output

**1b (host removal):** `/work_meta/preprocessing/02_host_removal/<sample>/`
- `*_R1_nonhost.fastq.gz`, `*_R2_nonhost.fastq.gz`, `*_singleton_nonhost.fastq.gz`
- `*_combined_nonhost.fastq.gz` — all three files concatenated

## Software
- BBDuk (BBMap 38.96) — adapter trimming (k=23, ktrim=r, trimq=20, minlen=50)
- BBMerge (BBMap 38.96) — merges overlapping read pairs
- Minimap2 2.28 — maps paired reads to host genome (`-ax sr`)
- Samtools 1.15.1 — filters unmapped reads (`-F 256 -f 4`) and converts to FASTQ

## Usage
```bash
# Step 1a: run serially over all samples
sbatch 1a_preprocess_qc.sh

# Step 1b: run as array (update N to match sample count)
sbatch --array=0-112 1b_host_removal.sh
```

## Notes
- Host reference: *Bombus impatiens* BIMP_2.2 (GCA_000188095.4); downloaded automatically if absent
- `1b` uses a file lock to prevent race conditions when building the Minimap2 index in parallel
- pigz is used for faster compression if available
