#!/bin/bash
#SBATCH --job-name=kraken2_screen
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/kraken2_screen_%A_%a.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/kraken2_screen_%A_%a.err

# ============================================================
# 08b_kraken2_screen.sh
# Author: Nickole Villabona
# Composition screen of the "non-host" reads: bacterial vs eukaryotic
# ============================================================
# Tests whether the rufocinctus contrast means what it appears to mean.
#
# The problem: host removal maps against *B. impatiens* BIMP_2.2 only, for
# all 22 Bombus species. Divergent lineages filter less host, so their
# "non-host" read counts are inflated by bee DNA that escaped the filter.
# rufocinctus is Cullumanobombus, one of the lineages furthest from
# impatiens, and it has 5.7M/4.7M non-host reads yet zero bins. If a large
# share of those reads are actually eukaryotic, then "depth is not limiting
# for rufocinctus" is wrong, multi-sample binning will not fix it, and the
# real fix is a closer host reference or a taxonomic pre-filter.
#
# Four samples, two per group:
#   rufocinctus (Cullumanobombus, far from impatiens, high read count)
#   Pyrobombus  (closer to impatiens, low read count) as the comparison
#
# Input:  ${READS_DIR}/<s>/<s>_R1/R2_nonhost.fastq.gz (read-only)
# Output: ${RESULTS}/kraken2/<s>.report  — full Kraken2 report
#         ${RESULTS}/kraken2/summary.tsv — appended domain-level fractions
#
# Tools: kraken2/2.1.2 (cluster module) + the capped standard-8 index
#        installed by 08a_kraken2_setup.sh
#
# Usage: sbatch 08a_kraken2_setup.sh          # once
#        sbatch --array=1-4 08b_kraken2_screen.sh
# ============================================================

set -euo pipefail

# Under sbatch the script runs from SLURM's spool directory, so BASH_SOURCE does
# NOT point at this directory (first run, 2026-11-25: all 22 tasks died with
# "config.sh: No such file or directory"). Use the absolute path; override with
# SCRIPT_DIR=... if the repo ever moves.
SCRIPT_DIR="${SCRIPT_DIR:-/dfs10/hammert-lab/nvillabo/Ch2/SYMBUS/scripts/12_multisample_binning}"
source "${SCRIPT_DIR}/config.sh"

module purge
module load kraken2/2.1.2

KRAKEN_DB=/dfs10/hammert-lab/nvillabo/databases/kraken2_standard8
OUT="${RESULTS}/kraken2"
SUMMARY="${OUT}/summary.tsv"

for f in hash.k2d opts.k2d taxo.k2d; do
    [[ -s "${KRAKEN_DB}/${f}" ]] || die "Kraken2 database incomplete (${KRAKEN_DB}/${f}) — run 08a_kraken2_setup.sh first"
done

# Two rufocinctus (depth_contrast) and two Pyrobombus with contrasting depth.
SCREEN_SAMPLES=(202008090837 202008090840 202007231044 202007251344)

: "${SLURM_ARRAY_TASK_ID:?Run with --array=1-4}"
SID="${SCREEN_SAMPLES[$(( SLURM_ARRAY_TASK_ID - 1 ))]:-}"
[[ -n "${SID}" ]] || { echo "[$(date)] No sample at index ${SLURM_ARRAY_TASK_ID} — exiting."; exit 0; }

R1="${READS_DIR}/${SID}/${SID}_R1_nonhost.fastq.gz"
R2="${READS_DIR}/${SID}/${SID}_R2_nonhost.fastq.gz"
[[ -s "${R1}" && -s "${R2}" ]] || die "missing non-host reads for ${SID}"

mkdir -p "${OUT}" "${LOG_DIR}"

REPORT="${OUT}/${SID}.report"
if [[ ! -s "${REPORT}" ]]; then
    echo "[$(date)] kraken2 on ${SID}"
    kraken2 \
      --db "${KRAKEN_DB}" \
      --threads "${SLURM_CPUS_PER_TASK}" \
      --paired --gzip-compressed \
      --report "${REPORT}.tmp" \
      --output /dev/null \
      "${R1}" "${R2}"
    mv "${REPORT}.tmp" "${REPORT}"
else
    echo "[$(date)] report exists for ${SID}, skipping kraken2"
fi

# ======== EXTRACT DOMAIN-LEVEL FRACTIONS ========
# Kraken2 report columns: pct, reads_clade, reads_taxon, rank, taxid, name
get_pct() {  # get_pct <rank-code> <name>
    awk -F'\t' -v r="$1" -v n="$2" '
      { gsub(/^ +/, "", $6) }
      $4 == r && $6 == n { printf "%.2f", $1; found=1; exit }
      END { if (!found) printf "0.00" }
    ' "${REPORT}"
}

UNCL=$(get_pct U unclassified)
BACT=$(get_pct D Bacteria)
EUK=$(get_pct D Eukaryota)
ARCH=$(get_pct D Archaea)
VIR=$(get_pct D Viruses)
# Residual host signal specifically: insects within Eukaryota.
ARTH=$(get_pct P Arthropoda)

echo
echo "[$(date)] ${SID} composition of the NON-HOST reads:"
printf '  unclassified %6s%%\n  Bacteria     %6s%%\n  Eukaryota    %6s%%  (Arthropoda %s%%)\n  Archaea      %6s%%\n  Viruses      %6s%%\n' \
  "${UNCL}" "${BACT}" "${EUK}" "${ARTH}" "${ARCH}" "${VIR}"

# Append to a shared summary. flock keeps concurrent array tasks from
# interleaving partial lines.
(
  flock -x 9
  [[ -s "${SUMMARY}" ]] || printf 'sample\tunclassified_pct\tbacteria_pct\teukaryota_pct\tarthropoda_pct\tarchaea_pct\tviruses_pct\n' > "${SUMMARY}"
  grep -qP "^${SID}\t" "${SUMMARY}" || \
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${SID}" "${UNCL}" "${BACT}" "${EUK}" "${ARTH}" "${ARCH}" "${VIR}" >> "${SUMMARY}"
) 9>"${SUMMARY}.lock"

echo "[$(date)] DONE ${SID} -> ${SUMMARY}"
