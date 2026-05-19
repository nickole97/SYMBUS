#!/bin/bash
#SBATCH --job-name=kofamscan_rerun
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs/kofamscan/kofamscan_rerun_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs/kofamscan/kofamscan_rerun_%A_%a.err

# Rerun for 11 MAGs that timed out in job 52097823 (original 12h limit)
# Increased to 24h. Usage: sbatch --array=1-11 8b_kofamscan_rerun.sh

DBCAN_DIR=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/08_annotation/dbcan_all
OUT_BASE=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/08_annotation/kofamscan
DB_DIR=/dfs10/hammert-lab/nvillabo/databases/kofamscan
MAG_LIST=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/08_annotation/kofamscan_missing_list.txt
LOG_DIR=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs/kofamscan

mkdir -p "$OUT_BASE" "$LOG_DIR"

module purge
module load mamba/24.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /dfs10/hammert-lab/nvillabo/envs/kofamscan_env

ENV_PATH=/dfs10/hammert-lab/nvillabo/envs/kofamscan_env
RUBY="${ENV_PATH}/bin/ruby"
EXEC_ANNOTATION="${ENV_PATH}/bin/exec_annotation"

find "${ENV_PATH}/lib/ruby" -name '*.gemspec' > /dev/null 2>&1 || true

echo "[$(date)] Using environment: kofamscan_env"
for attempt in 1 2 3; do
    if "$RUBY" "$EXEC_ANNOTATION" --version 2>/dev/null; then
        break
    fi
    echo "[$(date)] Ruby startup failed (attempt ${attempt}/3) — retrying in 20s..."
    sleep 20
    find "${ENV_PATH}/lib/ruby" -name '*.gemspec' > /dev/null 2>&1 || true
    if [ "$attempt" -eq 3 ]; then
        echo "[$(date)] ERROR: exec_annotation not found after 3 attempts"
        exit 1
    fi
done

MAG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MAG_LIST")
INPUT_FAA="${DBCAN_DIR}/${MAG}/uniInput.faa"
MAG_OUT="${OUT_BASE}/${MAG}"
TMP_DIR="${MAG_OUT}/tmp"

echo "[$(date)] MAG: $MAG (task ${SLURM_ARRAY_TASK_ID})"

[[ -f "$INPUT_FAA" ]] || { echo "[$(date)] ERROR: $INPUT_FAA not found"; exit 1; }

mkdir -p "$MAG_OUT" "$TMP_DIR"

"$RUBY" "$EXEC_ANNOTATION" \
  --cpu "$SLURM_CPUS_PER_TASK" \
  --ko-list  "${DB_DIR}/ko_list" \
  --profile  "${DB_DIR}/profiles" \
  --tmp-dir  "$TMP_DIR" \
  -f mapper \
  -o "${MAG_OUT}/ko_mapper.txt" \
  "$INPUT_FAA"

rm -rf "$TMP_DIR"

echo "[$(date)] Done — $(wc -l < "${MAG_OUT}/ko_mapper.txt") KO assignments written for $MAG"
echo "[$(date)] ==== End of KofamScan for $MAG ===="
