#!/bin/bash
#SBATCH --job-name=kofamscan_setup
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --partition=free
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs/kofamscan_setup_%A.out

# ============================================================
# 8b_kofamscan_setup.sh
# Author: Nickole Villabona
# One-time setup: install kofamscan_env and download KEGG DB
# ============================================================
# Creates the kofamscan_env conda environment and downloads
# the KofamScan database (ko_list + HMM profiles, ~8 GB).
# Must run on a compute node (never on login node).
# Also generates the MAG list needed for the job array.
#
# Input:  /work_meta/08_annotation/dbcan_all/  (to list MAGs)
# Output: /dfs10/hammert-lab/nvillabo/envs/kofamscan_env
#         /dfs10/hammert-lab/nvillabo/databases/kofamscan/
#         /work_meta/08_annotation/kofamscan_mag_list.txt
#
# Usage: sbatch 8b_kofamscan_setup.sh
# After done: sbatch --array=1-123 8b_kofamscan.sh
# ============================================================

ENV_PATH=/dfs10/hammert-lab/nvillabo/envs/kofamscan_env
DB_DIR=/dfs10/hammert-lab/nvillabo/databases/kofamscan
DBCAN_DIR=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/08_annotation/dbcan_all
LIST=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/08_annotation/kofamscan_mag_list.txt

mkdir -p "$DB_DIR"
mkdir -p /dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs

# ======== LOAD ENVIRONMENT ========
module purge
module load mamba/24.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"

# ======== CREATE CONDA ENV ========
if conda env list | grep -q "$ENV_PATH"; then
  echo "[$(date)] kofamscan_env already exists — skipping creation"
else
  echo "[$(date)] Creating kofamscan_env..."
  mamba create -y -p "$ENV_PATH" -c bioconda -c conda-forge kofamscan
  conda clean -a -f -y
  echo "[$(date)] kofamscan_env created at $ENV_PATH"
fi

conda activate "$ENV_PATH"
echo "[$(date)] Using environment: $ENV_PATH"
exec_annotation --version || { echo "[$(date)] ERROR: exec_annotation not found"; exit 1; }

# ======== DOWNLOAD DATABASE ========
if [[ -f "${DB_DIR}/ko_list" && -d "${DB_DIR}/profiles" ]]; then
  echo "[$(date)] KofamScan database already present — skipping download"
else
  echo "[$(date)] Downloading ko_list..."
  wget -q --show-progress -O "${DB_DIR}/ko_list.gz" https://www.genome.jp/ftp/db/kofam/ko_list.gz
  gunzip -f "${DB_DIR}/ko_list.gz"
  echo "[$(date)] ko_list ready"

  echo "[$(date)] Downloading profiles (~8 GB, this will take a while)..."
  wget -q --show-progress -O "${DB_DIR}/profiles.tar.gz" https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
  echo "[$(date)] Extracting profiles..."
  tar -xzf "${DB_DIR}/profiles.tar.gz" -C "$DB_DIR"
  rm "${DB_DIR}/profiles.tar.gz"
  echo "[$(date)] Database ready at $DB_DIR"
fi

# ======== GENERATE MAG LIST ========
ls "$DBCAN_DIR" | sort > "$LIST"
N=$(wc -l < "$LIST")
echo "[$(date)] MAG list written: $LIST ($N MAGs)"

echo "[$(date)] ==== Setup complete ===="
echo "[$(date)] Submit the job array with:"
echo "  sbatch --array=1-${N} 8b_kofamscan.sh"
