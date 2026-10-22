# Step 09 — Phylogenetic Analysis

Performs ANI-based species delineation for *Snodgrassella* MAGs (novel vs. known) and builds a multi-genome maximum-likelihood phylogeny using single-copy core orthologs (snod_phylo_v2).

## Inputs
**FastANI:**
- `/work_meta/06_MAGs_selected/*.fa` — focal Snodgrassella MAGs (lists hard-coded in script)

**snod_phylo_v2 pipeline:**
- Reference genomes downloaded from NCBI (see `snod_phylo_v2/scripts/00_download_refs.sh`)
- Our MAGs from `06_MAGs_selected/` (linked by `00b_link_our_MAGs.sh`)

## Outputs
**FastANI:**
- `/work_meta/09_phylogeny/fastani_snodgrassella/snod_allvsall.tsv` — all-vs-all ANI matrix

**snod_phylo_v2:**
- `03_prodigal/` — protein predictions per genome
- `04_orthofinder/` — single-copy core ortholog groups
- `06_alignments/` — MAFFT alignments, trimAl trimmed, concatenated supermatrix
- `07_iqtree/` — ML tree (`-m TEST`, 1000 bootstrap replicates)

## Software
- FastANI — pairwise ANI with `--matrix` output
- Prodigal — gene prediction (meta mode)
- OrthoFinder — single-copy core ortholog identification
- MAFFT — multiple sequence alignment
- trimAl — alignment trimming
- IQ-TREE — ML phylogeny (`-m TEST`)
- CondaEnv: `gtdbtk-2.5.2` (provides FastANI)

## Usage
```bash
# FastANI species delineation
sbatch snodgrassella_fastani.sh

# snod_phylo_v2 multi-step pipeline (run scripts in order)
bash snod_phylo_v2/scripts/00_download_refs.sh
bash snod_phylo_v2/scripts/00b_link_our_MAGs.sh
sbatch snod_phylo_v2/scripts/01_prodigal.sh
sbatch snod_phylo_v2/scripts/02_orthofinder.sh
python snod_phylo_v2/scripts/03_filter_core_OGs.py
sbatch snod_phylo_v2/scripts/04_align_concat.sh
sbatch snod_phylo_v2/scripts/05_iqtree.sh
```

## Notes
- FastANI script compares 17 novel Colombian *Snodgrassella* MAGs against 3 known *S. alvi* MAGs
- snod_phylo_v2 expands v1 to include Apis references, stingless bee genomes, and named species type strains (56 genomes total); see `snod_phylo_v2/README.md` for the full genome list
