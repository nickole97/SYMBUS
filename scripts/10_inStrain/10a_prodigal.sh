#!/bin/bash
#SBATCH --job-name=prodigal_instrain
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/logs/10a_prodigal_%j.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/logs/10a_prodigal_%j.err

###############################################################################
# 10a_prodigal.sh
# Gene calling con Prodigal sobre el CONCAT_FASTA renombrado
# IMPORTANTE: debe correr DESPUES de 10b_prepare_refs.sh
# Input:  10_inStrain/dereplicated_genomes.fasta (scaffolds renombrados)
# Output: 10_inStrain/dereplicated_genomes.genes.fna
#         10_inStrain/dereplicated_genomes.genes.faa
###############################################################################

# Cargar conda
module load miniconda3/23.5.2
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate binning_env

# Archivos
CONCAT_FASTA="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/dereplicated_genomes.fasta"
GENES_FNA="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/dereplicated_genomes.genes.fna"
GENES_FAA="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/dereplicated_genomes.genes.faa"

echo "=========================================="
echo "Starting Prodigal gene calling"
echo "Date: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "=========================================="

# Verificar que el CONCAT_FASTA existe
if [ ! -f "${CONCAT_FASTA}" ]; then
    echo "ERROR: CONCAT_FASTA no encontrado: ${CONCAT_FASTA}"
    echo "Corre primero 10b_prepare_refs.sh"
    exit 1
fi

echo "Input: ${CONCAT_FASTA}"
echo "Scaffolds totales: $(grep -c "^>" ${CONCAT_FASTA})"

# Correr Prodigal en modo meta sobre el CONCAT_FASTA completo
echo "[$(date)] Corriendo Prodigal en modo meta..."
prodigal \
    -i "${CONCAT_FASTA}" \
    -d "${GENES_FNA}" \
    -a "${GENES_FAA}" \
    -p meta \
    -q || true

echo ""
echo "=========================================="
echo "Prodigal completado!"
echo "Date: $(date)"
echo "Genes (fna): ${GENES_FNA}"
echo "Proteinas (faa): ${GENES_FAA}"
echo "Total genes predichos: $(grep -c "^>" ${GENES_FNA})"
echo "=========================================="
