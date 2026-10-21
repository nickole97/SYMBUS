#!/bin/bash
#SBATCH --job-name=dbcan_annotation
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --array=1-124
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/logs/dbcan_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/logs/dbcan_%A_%a.err

set -euo pipefail

echo "[$(date)] ==== Iniciando dbCAN SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} ===="

# === Rutas ===
MAG_DIR=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected_SD
DB_DIR=/dfs10/hammert-lab/nvillabo/databases/dbcan_db
OUT_DIR=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/08_annotation/dbcan_all
mkdir -p ${OUT_DIR}

# === Activar ambiente ===
module load miniconda3/23.5.2
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate dbcan5_env

# === Obtener MAG correspondiente a este array task ===
MAG=$(ls ${MAG_DIR}/*.fa | sed -n "${SLURM_ARRAY_TASK_ID}p")
SAMPLE=$(basename ${MAG} .fa)

echo "[$(date)] Procesando: ${SAMPLE}"

# === Correr dbCAN ===
run_dbcan CAZyme_annotation \
    --input_raw_data ${MAG} \
    --mode prok \
    --db_dir ${DB_DIR} \
    --output_dir ${OUT_DIR}/${SAMPLE} \
    --threads 8

echo "[$(date)] ==== ${SAMPLE} completado ===="
