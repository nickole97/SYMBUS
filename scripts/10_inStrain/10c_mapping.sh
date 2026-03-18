#!/bin/bash
#SBATCH --job-name=mapping_instrain
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --array=1-113
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/logs/10c_mapping_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/logs/10c_mapping_%A_%a.err

###############################################################################
# 10c_mapping.sh
# Mapeo de reads limpios contra MAGs dereplicated al 95% ANI
# Usando minimap2 + samtools (SLURM array, una muestra por job)
# Input:  preprocessing/clean_reads/*_R1_nonhost.fastq.gz
#         10_inStrain/dereplicated_genomes.fasta
# Output: 10_inStrain/10b_mapping/*.sorted.bam
###############################################################################

# Cargar conda
module load miniconda3/23.5.2
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate instrain_env

# Directorios
READS_DIR="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/preprocessing/clean_reads"
CONCAT_FASTA="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/dereplicated_genomes.fasta"
MAPPING_OUT="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/10b_mapping"

mkdir -p "${MAPPING_OUT}"

# Obtener lista de R1 reads y seleccionar el de este job
mapfile -t FILES < <(ls -1 "${READS_DIR}"/*_R1_nonhost.fastq.gz | sort)
NFILES=${#FILES[@]}
IDX=$((SLURM_ARRAY_TASK_ID-1))

if [[ $IDX -lt 0 || $IDX -ge $NFILES ]]; then
    echo "Task ID fuera de rango: ${SLURM_ARRAY_TASK_ID} (NFILES=${NFILES})"
    exit 0
fi

# Definir R1, R2 y sample ID
R1=${FILES[$IDX]}
R2=$(echo "${R1}" | sed 's/_R1_nonhost/_R2_nonhost/')
SAMPLE=$(basename "${R1}" _R1_nonhost.fastq.gz)

echo "=========================================="
echo "Mapping sample: ${SAMPLE}"
echo "Date: $(date)"
echo "Job ID: ${SLURM_JOB_ID}, Array task: ${SLURM_ARRAY_TASK_ID}"
echo "=========================================="
echo "R1: ${R1}"
echo "R2: ${R2}"

# Verificar que existen los reads
if [ ! -f "${R1}" ]; then
    echo "ERROR: R1 no encontrado: ${R1}"
    exit 1
fi

if [ ! -f "${R2}" ]; then
    echo "ERROR: R2 no encontrado: ${R2}"
    exit 1
fi

BAM_FILE="${MAPPING_OUT}/${SAMPLE}.sorted.bam"

# Saltar si ya existe
if [ -f "${BAM_FILE}" ] && [ -f "${BAM_FILE}.bai" ]; then
    echo "BAM ya existe para ${SAMPLE}, saltando..."
    exit 0
fi

# Mapear con minimap2 y ordenar con samtools
echo "[$(date)] Mapeando ${SAMPLE}..."
minimap2 \
    -ax sr \
    -t "${SLURM_CPUS_PER_TASK}" \
    "${CONCAT_FASTA}" \
    "${R1}" "${R2}" | \
samtools view \
    -@ "${SLURM_CPUS_PER_TASK}" \
    -bS - | \
samtools sort \
    -@ "${SLURM_CPUS_PER_TASK}" \
    -m 2G \
    -o "${BAM_FILE}" -

# Indexar BAM
echo "[$(date)] Indexando BAM..."
samtools index -@ "${SLURM_CPUS_PER_TASK}" "${BAM_FILE}"

# Estadisticas de mapeo
echo "[$(date)] Generando estadisticas..."
samtools flagstat "${BAM_FILE}" > "${BAM_FILE}.flagstat"

echo ""
echo "=========================================="
echo "Mapping completado para: ${SAMPLE}"
echo "Date: $(date)"
echo "BAM: ${BAM_FILE}"
cat "${BAM_FILE}.flagstat" | grep "mapped ("
echo "=========================================="
