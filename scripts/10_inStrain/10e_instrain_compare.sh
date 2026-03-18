#!/bin/bash
#SBATCH --job-name=instrain_compare
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/logs/10e_compare_%j.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/logs/10e_compare_%j.err

###############################################################################
# 10e_instrain_compare.sh
# Compara todos los perfiles inStrain entre si
# Calcula popANI entre todas las muestras por MAG
# Input:  10_inStrain/10c_profiles/*.IS/
#         10_inStrain/dereplicated_genomes.stb
# Output: 10_inStrain/10d_compare/
###############################################################################

# Cargar conda
module load miniconda3/23.5.2
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate instrain_env

# Directorios
PROFILE_OUT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/10c_profiles"
COMPARE_OUT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/10d_compare"
STB_FILE="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/dereplicated_genomes.stb"

mkdir -p "${COMPARE_OUT}"

echo "=========================================="
echo "inStrain compare - todas las muestras"
echo "Date: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "=========================================="

# Recopilar todos los perfiles .IS
mapfile -t PROFILES < <(ls -d "${PROFILE_OUT}"/*.IS)
N_PROFILES=${#PROFILES[@]}

echo "Perfiles encontrados: ${N_PROFILES}"

if [ ${N_PROFILES} -lt 2 ]; then
    echo "ERROR: Se necesitan al menos 2 perfiles para comparar"
    exit 1
fi

# Correr inStrain compare
echo "[$(date)] Corriendo inStrain compare..."
inStrain compare \
    -i "${PROFILES[@]}" \
    -o "${COMPARE_OUT}/all_samples" \
    -p "${SLURM_CPUS_PER_TASK}" \
    -s "${STB_FILE}" \
    --database_mode || echo "Compare terminó con advertencias"

echo ""
echo "=========================================="
echo "inStrain compare completado!"
echo "Date: $(date)"
echo "Output: ${COMPARE_OUT}"
ls "${COMPARE_OUT}/all_samples/output/" 2>/dev/null
echo "=========================================="
