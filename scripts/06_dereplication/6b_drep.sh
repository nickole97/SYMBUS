#!/bin/bash
#SBATCH --job-name=dRep99
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=18:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs/drep/drep99_%A.out

# ============================================================
# 6b_drep.sh
# Author: Nickole Villabona
# Dereplicate MAGs at 99% ANI (species-level)
# ============================================================
# Dereplicates selected MAGs to retain one representative genome
# per species (99% ANI threshold). Uses MASH for primary clustering
# and fastANI for refined secondary clustering.
#
# Input:  /work_meta/06_MAGs_selected/*.fa
# Output: /work_meta/07_dRep_99_v3/
#         - dereplicated_genomes/*.fa (one representative per species)
#         - data_tables/ (clustering and quality information)
#         - figures/ (dendrograms and plots)
#
# Tools: dRep (primary ANI 90%, secondary ANI 99%, single-linkage)
#
# Usage: sbatch 6b_drep.sh
# Runtime: ~6-12 hours depending on number of MAGs
# ============================================================

# ======== MODULES ========
module load mamba/24.3.0
conda activate drep

# ======== PATHS ========
MAG_DIR="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected"
OUTPUT_DIR="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/07_dRep_99_v3"

echo "=========================================="
echo "Starting dRep dereplication at 99% ANI"
echo "=========================================="
echo "Date: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "MAG directory: ${MAG_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Number of MAGs: $(ls ${MAG_DIR}/*.fa | wc -l)"
echo "=========================================="
echo ""

# Create output directory
mkdir -p ${OUTPUT_DIR}

# ======== RUN DREP ========
dRep dereplicate ${OUTPUT_DIR} \
  -g ${MAG_DIR}/*.fa \
  -pa 0.9 \              # Primary ANI: 90% (MASH passes genomes with ≥90% similarity)
  -sa 0.99 \             # Secondary ANI: 99% (fastANI clusters genomes at ≥99% ANI)
  -nc 16 \               # Number of CPUs
  --S_algorithm fastANI \  # Force fastANI for secondary clustering
  -clusterAlg single     # Single-linkage: avoid unnecessary cluster subdivision

echo ""
echo "=========================================="
echo "dRep completed!"
echo "Date: $(date)"
echo "Number of dereplicated MAGs: $(ls ${OUTPUT_DIR}/dereplicated_genomes/*.fa 2>/dev/null | wc -l)"
echo "=========================================="
