#!/bin/bash
#SBATCH --job-name=minpath
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs/minpath/minpath_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs/minpath/minpath_%A_%a.err


# ============================================================
# 8d_minpath.sh
# Author: Nickole Villabona
# KEGG pathway reconstruction with MinPath
# ============================================================
# Infers which KEGG pathways are genuinely present in each MAG
# from its KofamScan KO assignments, using MinPath's parsimony
# approach (smallest set of pathways explaining the observed KOs).
#
# Input:  /work_meta/08_annotation/kofamscan/<MAG>/ko_mapper.txt
#         /work_meta/08_annotation/minpath_mag_list.txt (121 MAGs)
# Output: /work_meta/08_annotation/minpath/<MAG>/
#         - minpath.report   (one line per pathway; "minpath 1" = present)
#         - minpath.details  (KOs supporting each pathway)
#
# Tools: MinPath (github.com/mgtools/MinPath @ 46d3e81) + GLPK 4.6 solver
#        Installed manually at /dfs10/hammert-lab/nvillabo/tools/MinPath
#        (see README — no setup script; bioconda's minpath did not work)
#
# Usage: sbatch --array=1-121 8d_minpath.sh
# ============================================================
KOFAM_DIR=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/08_annotation/kofamscan
OUT_DIR=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/08_annotation/minpath
MAG_LIST=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/08_annotation/minpath_mag_list.txt
MINPATH=/dfs10/hammert-lab/nvillabo/tools/MinPath/MinPath.py

MAG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MAG_LIST")
KO_IN="${KOFAM_DIR}/${MAG}/ko_mapper.txt"
MAG_OUT="${OUT_DIR}/${MAG}"

echo "[$(date)] MAG: $MAG (task ${SLURM_ARRAY_TASK_ID})"

N_KOS=$(awk -F'\t' 'NF==2' "$KO_IN" | wc -l)
echo "[$(date)] KOs with assignment: $N_KOS"

if [ "$N_KOS" -eq 0 ]; then
    echo "[$(date)] WARNING: no KOs for $MAG — skipping"
    exit 0
fi

mkdir -p "$MAG_OUT"

# MinPath writes temp files (test.mps, test.mps.LPout) to cwd
cd "$MAG_OUT"

python3 "$MINPATH" \
    -ko      "$KO_IN" \
    -report  "${MAG_OUT}/minpath.report" \
    -details "${MAG_OUT}/minpath.details"

N_ACTIVE=$(grep -c 'minpath 1' "${MAG_OUT}/minpath.report" 2>/dev/null || echo 0)
echo "[$(date)] Done: $N_ACTIVE pathways active for $MAG"
