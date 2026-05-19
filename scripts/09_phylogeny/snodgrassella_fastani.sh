#!/bin/bash
#SBATCH --job-name=fastani_snod
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --partition=free
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/logs/fastani_snod_%A.out

# ============================================================
# FastANI — Snodgrassella species delineation
# Novel (Colombian, no-species assignment) vs. known species
# ============================================================

MAG_DIR=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected
OUT_DIR=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/09_phylogeny/fastani_snodgrassella
mkdir -p "$OUT_DIR"

module purge
module load mamba/24.3.0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /dfs10/hammert-lab/nvillabo/envs/gtdbtk-2.5.2

echo "[$(date)] FastANI version: $(fastANI --version 2>&1)"

# ---- Novel Snodgrassella (Colombian, no species in GTDB) ----
cat > "${OUT_DIR}/novel_snod.txt" << 'FILELIST'
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV006__54_sub.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV009__38.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV015__NV015_metabat.3.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV015__NV015_metabat.8.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV016__43.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV018__21_sub.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV020__NV020_maxbin.001.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV022__NV022_metabat.5.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV029__1.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV032__17_sub.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV033__14.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV034__NV034_maxbin.003.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV034__NV034_maxbin.005.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV048__NV048_metabat.2.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV057__NV057_metabat.2.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV059__NV059_metabat.2.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV059__NV059_metabat.4.fa
FILELIST

# ---- Known Snodgrassella from the same collection ----
cat > "${OUT_DIR}/known_snod.txt" << 'FILELIST'
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/202006130802__202006130802_maxbin.001.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/202006140655__39_sub.fa
/dfs10/hammert-lab/nvillabo/Ch2/work_meta/06_MAGs_selected/NV002__NV002_metabat.1.fa
FILELIST

# ---- All Snodgrassella combined ----
cat "${OUT_DIR}/novel_snod.txt" "${OUT_DIR}/known_snod.txt" > "${OUT_DIR}/all_snod.txt"

echo "[$(date)] Running all-vs-all FastANI on all Snodgrassella MAGs..."
fastANI \
    --ql "${OUT_DIR}/all_snod.txt" \
    --rl "${OUT_DIR}/all_snod.txt" \
    -o "${OUT_DIR}/snod_allvsall.tsv" \
    -t "${SLURM_CPUS_PER_TASK}" \
    --matrix

echo "[$(date)] Done. Output: ${OUT_DIR}/snod_allvsall.tsv"
echo "[$(date)] Results preview:"
sort -k3 -rn "${OUT_DIR}/snod_allvsall.tsv" | head -20

echo "[$(date)] ==== FastANI Snodgrassella complete ===="
