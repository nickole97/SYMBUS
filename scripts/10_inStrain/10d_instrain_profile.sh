#!/bin/bash
#SBATCH --job-name=instrain_profile
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --array=1-113
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/logs/10d_profile_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/logs/10d_profile_%A_%a.err

###############################################################################
# 10d_instrain_profile.sh
# Corre inStrain profile para cada muestra
# Input:  10_inStrain/10b_mapping/*.sorted.bam
#         10_inStrain/dereplicated_genomes.fasta
#         10_inStrain/dereplicated_genomes.genes.fna
#         10_inStrain/dereplicated_genomes.stb
# Output: 10_inStrain/10c_profiles/*.IS/
###############################################################################

# Cargar conda
module load miniconda3/23.5.2
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate instrain_env

# Directorios y archivos de referencia
MAPPING_OUT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/10b_mapping"
PROFILE_OUT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/10c_profiles"
CONCAT_FASTA="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/dereplicated_genomes.fasta"
GENES_FNA="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/dereplicated_genomes.genes.fna"
STB_FILE="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/dereplicated_genomes.stb"

mkdir -p "${PROFILE_OUT}"

# Obtener lista de BAMs y seleccionar el de este job
mapfile -t FILES < <(ls -1 "${MAPPING_OUT}"/*.sorted.bam | sort)
NFILES=${#FILES[@]}
IDX=$((SLURM_ARRAY_TASK_ID-1))

if [[ $IDX -lt 0 || $IDX -ge $NFILES ]]; then
    echo "Task ID fuera de rango: ${SLURM_ARRAY_TASK_ID} (NFILES=${NFILES})"
    exit 0
fi

BAM=${FILES[$IDX]}
SAMPLE=$(basename "${BAM}" .sorted.bam)
PROFILE_DIR="${PROFILE_OUT}/${SAMPLE}.IS"

echo "=========================================="
echo "inStrain profile: ${SAMPLE}"
echo "Date: $(date)"
echo "Job ID: ${SLURM_JOB_ID}, Array task: ${SLURM_ARRAY_TASK_ID}"
echo "=========================================="

# Saltar si ya existe
if [ -d "${PROFILE_DIR}" ]; then
    echo "Perfil ya existe para ${SAMPLE}, saltando..."
    exit 0
fi

# Correr inStrain profile
echo "[$(date)] Corriendo inStrain profile..."
inStrain profile \
    "${BAM}" \
    "${CONCAT_FASTA}" \
    -o "${PROFILE_DIR}" \
    -p "${SLURM_CPUS_PER_TASK}" \
    -g "${GENES_FNA}" \
    -s "${STB_FILE}" \
    -l 0.95 \
    -c 1 \
    -f 0.05 \

echo ""
echo "=========================================="
echo "inStrain profile completado: ${SAMPLE}"
echo "Date: $(date)"
echo "Output: ${PROFILE_DIR}"
echo "=========================================="
