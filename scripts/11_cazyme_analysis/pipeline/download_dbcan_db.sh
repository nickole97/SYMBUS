#!/bin/bash
#SBATCH --job-name=dbcan_download
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/logs/dbcan_download_%j.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/logs/dbcan_download_%j.err

set -euo pipefail

echo "[$(date)] ==== Iniciando descarga de bases de datos dbCAN2 ===="

module load miniconda3/23.5.2
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dbcan5_env

DB_DIR=/dfs10/hammert-lab/nvillabo/databases/dbcan_db
mkdir -p ${DB_DIR}

run_dbcan database --db_dir ${DB_DIR} --aws_s3

echo "[$(date)] ==== Todo listo! ===="