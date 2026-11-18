#!/bin/bash
#SBATCH --job-name=ms_depth
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/depth_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/depth_%A_%a.err

# ============================================================
# 02_depth_matrix.sh
# Author: Nickole Villabona
# Build the multi-sample depth products for one assembly
# ============================================================
# One array task per assembly (1..N_SUBSET). Consumes every cross-mapped
# BAM for that assembly and produces the three different coverage formats
# the three binners need:
#
#   depth.txt          multi-column jgi table      -> MetaBAT2 (native)
#   abund_<sample>.txt 2-column (contig, depth)    -> MaxBin2 via -abund_list
#   coverage_table.tsv per 10kb chunk, all samples -> CONCOCT
#
# This is deliberately the ONLY step that reads the BAMs, so that deleting
# them afterwards is safe. Do not move the CONCOCT coverage table into
# 03c: it would need the BAMs back after they are gone.
#
# NOTE on MaxBin2: the single-sample pipeline passes the jgi depth file to
# MaxBin2's -abund, which expects two columns (contig, abundance) and so
# reads contigLen as abundance. That is bug (a) in the project bug ledger
# and is NOT fixed in the existing results. Here MaxBin2 gets correctly
# formatted 2-column files from the start.
#
# Input:  ${SCRATCH}/01_crossmap/<asm>/*.sorted.bam
#         ${ASSEMBLY_DIR}/<asm>/final.contigs.fa (read-only)
# Output: ${SCRATCH}/02_depth/<asm>/
#         - depth.txt, bam_order.txt
#         - abund_<sample>.txt, abund_list.txt
#         - contigs_10K.fa, contigs_10K.bed, coverage_table.tsv
#
# Tools: jgi_summarize_bam_contig_depths (quast_env — same binary that made
#        the original single-sample depth files; see config.sh)
#        cut_up_fasta.py + concoct_coverage_table.py (binning_env)
#        Samtools 1.15.1 (module; concoct_coverage_table.py shells out to it)
#
# Usage: sbatch --array=1-22%20 02_depth_matrix.sh
# ============================================================

set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"
check_subset

module purge
module load samtools/1.15.1

THREADS="${SLURM_CPUS_PER_TASK:-8}"

: "${SLURM_ARRAY_TASK_ID:?Run with --array=1-${N_SUBSET}}"
ASM="$(sample_from_index "${SLURM_ARRAY_TASK_ID}")"
[[ -n "${ASM}" ]] || { echo "[$(date)] No sample at index ${SLURM_ARRAY_TASK_ID} — exiting."; exit 0; }

CONTIGS="${ASSEMBLY_DIR}/${ASM}/final.contigs.fa"
BAM_DIR="${SCRATCH}/01_crossmap/${ASM}"
OUT="${SCRATCH}/02_depth/${ASM}"
mkdir -p "${OUT}" "${LOG_DIR}"

echo "[$(date)] assembly ${ASM} (task ${SLURM_ARRAY_TASK_ID}/${N_SUBSET})"

# ======== GATHER BAMs IN A STABLE ORDER ========
# Order follows backmap_pairs.tsv, not shell globbing, so the column order of
# depth.txt is reproducible and matches bam_order.txt.
BAMS=()
while read -r RD; do
    B="${BAM_DIR}/${ASM}__${RD}.sorted.bam"
    [[ -s "${B}" ]] || die "missing BAM: ${B} — rerun the failed 01_crossmap tasks"
    [[ -s "${B}.bai" ]] || die "missing BAM index: ${B}.bai"
    BAMS+=("${B}")
done < <(awk -F'\t' -v a="${ASM}" '$1==a {print $2}' "${PAIRS_FILE}")

N_BAMS="${#BAMS[@]}"
[[ "${N_BAMS}" -gt 0 ]] || die "no BAMs found for ${ASM}"
printf '%s\n' "${BAMS[@]}" > "${OUT}/bam_order.txt"
echo "[$(date)] ${N_BAMS} BAMs for ${ASM}"

# ======== 1) MULTI-COLUMN DEPTH (MetaBAT2) ========
DEPTH="${OUT}/depth.txt"
if [[ ! -s "${DEPTH}" ]]; then
    echo "[$(date)] jgi_summarize_bam_contig_depths over ${N_BAMS} BAMs"
    "${JGI_DEPTH}" --outputDepth "${DEPTH}.tmp" "${BAMS[@]}"
    mv "${DEPTH}.tmp" "${DEPTH}"
else
    echo "[$(date)] depth.txt exists, skipping"
fi

# jgi layout: contigName, contigLen, totalAvgDepth, then per BAM a mean and a
# variance column. Verify before deriving anything from column positions.
EXPECTED_COLS=$(( 3 + 2 * N_BAMS ))
ACTUAL_COLS=$(head -1 "${DEPTH}" | awk -F'\t' '{print NF}')
[[ "${ACTUAL_COLS}" -eq "${EXPECTED_COLS}" ]] || \
    die "depth.txt has ${ACTUAL_COLS} columns, expected ${EXPECTED_COLS} for ${N_BAMS} BAMs"
N_CONTIGS=$(( $(wc -l < "${DEPTH}") - 1 ))
[[ "${N_CONTIGS}" -gt 0 ]] || die "depth.txt has no contig rows"
echo "[$(date)] depth.txt OK: ${N_CONTIGS} contigs x ${N_BAMS} samples"

# ======== 2) PER-SAMPLE 2-COLUMN ABUNDANCE (MaxBin2) ========
ABUND_LIST="${OUT}/abund_list.txt"
: > "${ABUND_LIST}"
i=0
while read -r RD; do
    i=$(( i + 1 ))
    COL=$(( 2 * i + 2 ))          # mean-depth column for BAM i
    AB="${OUT}/abund_${RD}.txt"
    awk -F'\t' -v c="${COL}" 'NR>1 {printf "%s\t%s\n", $1, $c}' "${DEPTH}" > "${AB}"
    [[ -s "${AB}" ]] || die "empty abundance file: ${AB}"
    echo "${AB}" >> "${ABUND_LIST}"
done < <(awk -F'\t' -v a="${ASM}" '$1==a {print $2}' "${PAIRS_FILE}")

[[ "$(wc -l < "${ABUND_LIST}")" -eq "${N_BAMS}" ]] || \
    die "abund_list.txt has $(wc -l < "${ABUND_LIST}") entries, expected ${N_BAMS}"
echo "[$(date)] ${N_BAMS} MaxBin2 abundance files written"

# ======== 3) CUT-UP CONTIGS + COVERAGE TABLE (CONCOCT) ========
CUT_FA="${OUT}/contigs_10K.fa"
CUT_BED="${OUT}/contigs_10K.bed"
COV="${OUT}/coverage_table.tsv"

if [[ ! -s "${CUT_FA}" || ! -s "${CUT_BED}" ]]; then
    echo "[$(date)] cutting contigs into ${CONCOCT_CHUNK} bp chunks"
    "${CONCOCT_CUTUP}" "${CONTIGS}" -c "${CONCOCT_CHUNK}" -o 0 --merge_last \
        -b "${CUT_BED}.tmp" > "${CUT_FA}.tmp"
    mv "${CUT_BED}.tmp" "${CUT_BED}"
    mv "${CUT_FA}.tmp"  "${CUT_FA}"
else
    echo "[$(date)] cut-up files exist, skipping"
fi

if [[ ! -s "${COV}" ]]; then
    echo "[$(date)] building CONCOCT coverage table over ${N_BAMS} BAMs"
    "${CONCOCT_COVTAB}" "${CUT_BED}" "${BAMS[@]}" > "${COV}.tmp"
    mv "${COV}.tmp" "${COV}"
else
    echo "[$(date)] coverage_table.tsv exists, skipping"
fi

# contig column + one column per sample
COV_COLS=$(head -1 "${COV}" | awk -F'\t' '{print NF}')
[[ "${COV_COLS}" -eq $(( N_BAMS + 1 )) ]] || \
    die "coverage_table.tsv has ${COV_COLS} columns, expected $(( N_BAMS + 1 ))"
[[ "$(wc -l < "${COV}")" -gt 1 ]] || die "coverage_table.tsv has no data rows"
echo "[$(date)] coverage_table.tsv OK: $(( $(wc -l < "${COV}") - 1 )) chunks x ${N_BAMS} samples"

# ======== 4) RELEASE THE BAMs ========
# Only now, with all three products validated above. Peak scratch stays at
# roughly one assembly's worth of BAMs instead of all of them.
if [[ "${DELETE_BAMS_AFTER_DEPTH}" -eq 1 ]]; then
    FREED=$(du -sh "${BAM_DIR}" 2>/dev/null | cut -f1)
    echo "[$(date)] deleting ${N_BAMS} BAMs for ${ASM} (${FREED}); DELETE_BAMS_AFTER_DEPTH=1"
    rm -rf "${BAM_DIR}"
    # The minimap2 index is no longer needed either once mapping is done.
    rm -f "${SCRATCH}/00_index/${ASM}.mmi"
else
    echo "[$(date)] keeping BAMs (DELETE_BAMS_AFTER_DEPTH=0)"
fi

echo "[$(date)] DONE ${ASM}"
