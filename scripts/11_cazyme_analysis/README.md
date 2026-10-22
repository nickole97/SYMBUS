# Step 11 — CAZyme Annotation & Analysis

Annotates CAZymes in all MAGs with dbCAN5, consolidates results into a MAG × CAZyme family count matrix, and performs comparative R-based analyses (diversity, heatmaps, ordination) stratified by bacterial genus, host species, and geography.

## Inputs
**dbCAN pipeline:**
- `/work_meta/06_MAGs_selected_SD/*.fa` — selected MAGs

**R analysis (`cazyme_analysis.Rmd`):**
- `cazyme_matrix.tsv` — MAG × CAZyme family matrix (high-quality MAGs only)
- `cazyme_matrix_all.tsv` — full matrix (all MAGs)
- `combined_taxonomy_splitted.txt` — GTDB taxonomy table with split lineage columns
- `genome_to_taxonomy.tsv` — MAG-to-taxonomy mapping
- `bombus_sample_metadata_2026.csv` — sample metadata (host species, country, collection date)

## Outputs
**dbCAN pipeline:**
- `/work_meta/08_annotation/dbcan_all/<MAG>/overview.tsv` — per-MAG CAZyme hits
- `/work_meta/08_annotation/dbcan_all/<MAG>/uniInput.faa` — Prodigal protein predictions (reused by KofamScan)
- Consolidated matrix produced by `consolidate_cazymes.py`

**R analysis:**
- HTML report with CAZyme overview plots, heatmaps, PCoA/NMDS ordinations, UpSet plots, and genus-level comparisons

## Software
**Pipeline:**
- dbCAN5 (`run_dbcan CAZyme_annotation --mode prok`) — HMMER + DIAMOND + dbCAN_sub annotation
- Python 3 / pandas (`consolidate_cazymes.py`) — filters hits with ≥ 2 tools and builds count matrix
- CondaEnv: `dbcan5_env`

**R packages:**
- tidyverse, ggplot2, patchwork — data wrangling and plotting
- pheatmap — clustered heatmaps
- vegan — Bray-Curtis dissimilarity, PCoA/NMDS ordination
- ggrepel — non-overlapping labels
- scales — axis formatting
- UpSetR — set intersection plots

## Usage
```bash
# One-time database download
sbatch pipeline/download_dbcan_db.sh

# dbCAN annotation array
sbatch pipeline/dbcan_all.sh          # --array=1-124 set inside script

# Consolidate into matrix
python pipeline/consolidate_cazymes.py \
  --input_dir /path/to/dbcan_all \
  --output cazyme_matrix.tsv

# R analysis (render locally or on HPC with R)
Rscript -e "rmarkdown::render('cazyme_analysis.Rmd')"
```

## Notes
- `consolidate_cazymes.py` defaults to `--min_tools 2` (high-confidence hits); use `--collapse_subfamilies` to merge subfamilies (e.g., GH43_4 → GH43)
- dbCAN database downloaded from AWS S3 (`--aws_s3` flag)
- The Rmd joins the CAZyme matrix with taxonomy and metadata by sample ID extracted from the MAG name prefix (everything before `__`)
