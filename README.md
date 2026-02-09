# SYMBUS

**Symbiosis Between Bumblebees and their gut microbiomes Using Metagenomics**

SYMBUS is a repository that contains all the scripts and workflows used to perform the analysis of the symbiosis between bumblebees and their gut microbiomes from the metagenomes of *Bombus* species.

## Overview

This pipeline processes metagenomic sequencing data from bumblebee gut samples (both US and Colombian populations) to reconstruct, annotate, and quantify microbial genomes (MAGs - Metagenome-Assembled Genomes). The workflow starts from raw read quality control through taxonomic and functional annotation. Some steps have already been implemented, other are currently under development and will be incorporated soon. 

## Pipeline Workflow

This diagram describes the complete workflow of the project:
<img src="https://github.com/user-attachments/assets/b53eed85-addf-4831-9bf0-5de6e3f891af" alt="SYMBUS pipeline" width="800">

The pipeline consists of eight major steps:
1. **Preprocessing**: Quality filtering and host genome removal
2. **Assembly**: Metagenomic assembly using MEGAHIT
3. **Mapping**: Read mapping back to assembled contigs
4. **Binning**: Genome binning using three complementary tools
5. **Quality Control**: MAG quality assessment with CheckM
6. **Dereplication**: Selection and dereplication of high-quality MAGs
7. **Abundance**: Quantification of MAG abundance across samples
8. **Annotation**: Taxonomic and functional annotation

## Computing Environment

All analyses were performed on the [HPC3 cluster](https://rcic.uci.edu/about/intro.html) at the **University of California, Irvine** running Linux OS (Ubuntu 24).

## Requirements

### Software Dependencies

The pipeline uses the following bioinformatics tools:

**Quality Control & Preprocessing:**
- [BBMap/BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) (v38.96) - adapter trimming and quality filtering
- [BBMerge](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) (v38.96) - read merging

**Host Removal & Mapping:**
- [Minimap2](https://github.com/lh3/minimap2) (v2.28) - fast sequence alignment
- [Samtools](http://www.htslib.org/) (v1.15.1) - BAM file manipulation

**Assembly:**
- [MEGAHIT](https://github.com/voutcn/megahit) (v1.2.9) - metagenomic assembler

**Binning:**
- [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/) - binning based on coverage and composition
- [MaxBin2](https://sourceforge.net/projects/maxbin2/) - binning using abundance and tetranucleotide frequency
- [CONCOCT](https://github.com/BinPro/CONCOCT) - binning using coverage and composition
- [DAS Tool](https://github.com/cmks/DAS_Tool) - bin refinement and dereplication across binners

**Quality Assessment:**
- [CheckM](https://github.com/Ecogenomics/CheckM) - assessment of genome completeness and contamination
- [dRep](https://github.com/MrOlm/drep) - genome dereplication at species level (99% ANI)

**Abundance Estimation:**
- [CoverM](https://github.com/wwood/CoverM) - genome coverage and abundance calculation

**Annotation:**
- [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) - taxonomic classification using GTDB database
- [Prokka](https://github.com/tseemann/prokka) - functional annotation (planned)
- [KOfamScan](https://github.com/takaram/kofam_scan) - KEGG pathway annotation (planned)
- [dbCAN2](https://github.com/linnabrown/run_dbcan) - CAZyme annotation (planned)

### Conda Environments

The pipeline uses several conda environments to manage dependencies:
- `bbmap/38.96` (module)
- `binning_env` (MetaBAT2, MaxBin2, CONCOCT, DAS Tool)
- `checkm_env` (CheckM)
- `drep` (dRep)
- `coverm_env` (CoverM)
- `quast_env` (MetaQUAST)

## Repository Structure

```
SYMBUS/
├── README.md                          # This file
├── assets/                            # Images and diagrams
│   └── SYMBUS_PIPELINE.png
├── scripts/                           # Analysis scripts
│   ├── 1_preprocessing/               # Quality control and host removal
│   ├── 2_assembly/                    # Metagenomic assembly
│   ├── 3_mapping/                     # Read mapping to contigs
│   ├── 4_binning/                     # Genome binning
│   ├── 5_quality_control/             # MAG quality assessment
│   ├── 6_dereplication/               # MAG selection and dereplication
│   ├── 7_abundance/                   # Abundance estimation
│   └── 8_annotation/                  # Taxonomic and functional annotation
└── notebooks/
    └── 00_pipeline_notebook.md        # Analysis notebook with all steps
```

## Usage

### General Workflow

Scripts are organized by major workflow steps (numbered 1-8) with substeps labeled alphabetically (a, b, c). **Scripts should be run sequentially** following the numbered order.

Each script is designed to run as an **SLURM array job** on the HPC cluster, processing multiple samples in parallel. The scripts use a common `samples.txt` file (generated during preprocessing) to iterate through all samples.

### Running the Pipeline

1. **Start with preprocessing** (`1_preprocessing/`) to quality-filter reads and remove host contamination
2. **Proceed through assembly** (`2_assembly/`) to generate contigs for each sample
3. **Continue sequentially** through mapping, binning, quality control, dereplication, and abundance estimation
4. **Finish with annotation** to assign taxonomy and function to your MAGs

Each script contains detailed header comments explaining:
- Purpose and inputs/outputs
- Required modules and dependencies  
- Key parameters
- Expected runtime
- Example SLURM submission commands

See individual scripts for specific usage and parameters.

### Key Quality Thresholds

MAGs are selected based on these criteria:
- **High quality**: Completeness ≥90%, Contamination <5%
- **Medium quality**: Completeness ≥70%, Contamination <10%
- **Strain diversity**: Coverage ≥5×, Breadth ≥50%
- **Dereplication**: 99% ANI threshold (species-level)

## Input Data

The pipeline expects paired-end Illumina metagenomic reads:
- US bumblebee samples: 91 samples from Western US *Bombus* species
- Colombian bumblebee samples: Additional samples from Colombian populations
- Format: Interleaved or paired FASTQ files (`.fq.gz` or `.fastq.gz`)

## Output Data

Key outputs include:
- Quality-filtered, host-removed reads
- Assembled contigs (minimum 1000 bp)
- Genome bins from individual samples
- High-quality, dereplicated MAGs (species representatives)
- Taxonomic assignments (GTDB taxonomy)
- Abundance matrices across all samples
- Functional annotations (in progress)

## Contributors

- **Nickole Villabona** - Pipeline development and analysis

## Contact

For questions or issues, please open an issue on this repository or contact:
- Nickole Villabona: [nvillabo@uci.edu]
---

