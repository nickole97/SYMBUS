#!/bin/bash
#SBATCH --job-name=prepare_refs
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/logs/10b_prepare_refs_%j.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/logs/10b_prepare_refs_%j.err

###############################################################################
# 10b_prepare_refs.sh
# Prepara archivos de referencia para inStrain:
#   1. Concatena MAGs dereplicated renombrando scaffolds (evita duplicados)
#   2. Crea el scaffold-to-bin (STB) file
# Input:  07_dRep_95_correct/dereplicated_genomes/*.fa
# Output: 10_inStrain/dereplicated_genomes.fasta
#         10_inStrain/dereplicated_genomes.stb
###############################################################################

# Directorios
DEREP_GENOMES="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/07_dRep_95_correct/dereplicated_genomes"
CONCAT_FASTA="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/dereplicated_genomes.fasta"
STB_FILE="/dfs10/hammert-lab/nvillabo/Ch2/work_meta/10_inStrain/dereplicated_genomes.stb"

echo "=========================================="
echo "Preparing reference files for inStrain"
echo "Date: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "=========================================="

# PASO 1: Concatenar MAGs renombrando scaffolds para evitar duplicados
# Formato nuevo: >BINNAME__scaffold_original
echo "Creando CONCAT_FASTA con scaffolds renombrados..."
> "${CONCAT_FASTA}"

for genome in "${DEREP_GENOMES}"/*.fa; do
    bin_name=$(basename "${genome}" .fa)
    # Renombrar cada header: >scaffold -> >BINNAME__scaffold
    awk -v bin="${bin_name}" '
        /^>/ { print ">" bin "__" substr($0, 2) }
        !/^>/ { print }
    ' "${genome}" >> "${CONCAT_FASTA}"
done

echo "CONCAT_FASTA creado: ${CONCAT_FASTA}"
echo "Total scaffolds: $(grep -c "^>" ${CONCAT_FASTA})"

# PASO 2: Crear STB file a partir del CONCAT_FASTA renombrado
# Formato STB: scaffold_name \t bin_name
echo ""
echo "Creando STB file..."
> "${STB_FILE}"

for genome in "${DEREP_GENOMES}"/*.fa; do
    bin_name=$(basename "${genome}" .fa)
    grep "^>" "${genome}" | sed 's/^>//' | awk -v bin="${bin_name}" '{
        print bin"__"$1"\t"bin
    }' >> "${STB_FILE}"
done

echo "STB file creado: ${STB_FILE}"
echo "Total scaffolds en STB: $(wc -l < ${STB_FILE})"

echo ""
echo "=========================================="
echo "Referencias listas!"
echo "Date: $(date)"
echo "  CONCAT_FASTA: ${CONCAT_FASTA}"
echo "  STB_FILE:     ${STB_FILE}"
echo "=========================================="
