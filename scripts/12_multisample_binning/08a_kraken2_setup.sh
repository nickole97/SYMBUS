#!/bin/bash
#SBATCH --job-name=kraken2_setup
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --partition=standard
#SBATCH --account=hammert_lab
#SBATCH --output=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/kraken2_setup_%j.out
#SBATCH --error=/dfs10/hammert-lab/nvillabo/Ch2/work_meta/12_multisample_test/logs/kraken2_setup_%j.err

# ============================================================
# 08a_kraken2_setup.sh
# Author: Nickole Villabona
# One-time: fetch a Kraken2 database for the read-composition screen
# ============================================================
# kraken2/2.1.2 is available as a cluster module, but no Kraken2 database
# exists under databases/ — this downloads one.
#
# Uses the CAPPED "standard-8" index (~8 GB) rather than the full standard
# index (~70 GB). The screen in 08b only needs domain-level resolution
# (bacterial vs eukaryotic/host), which the capped index gives at a fraction
# of the download and memory. Do NOT use this index for anything needing
# species-level read assignment.
#
# Output: ${KRAKEN_DB}/{hash.k2d,opts.k2d,taxo.k2d}
#
# Usage: sbatch 08a_kraken2_setup.sh
# Note: run once. 08b checks for the database and tells you to run this.
# ============================================================

set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

KRAKEN_DB=/dfs10/hammert-lab/nvillabo/databases/kraken2_standard8
DB_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz"
TARBALL="${KRAKEN_DB}/k2_standard_08gb.tar.gz"

mkdir -p "${KRAKEN_DB}" "${LOG_DIR}"

if [[ -s "${KRAKEN_DB}/hash.k2d" && -s "${KRAKEN_DB}/taxo.k2d" && -s "${KRAKEN_DB}/opts.k2d" ]]; then
    echo "[$(date)] database already present in ${KRAKEN_DB}:"
    du -sh "${KRAKEN_DB}"
    exit 0
fi

echo "[$(date)] downloading ${DB_URL}"
wget -q --show-progress -O "${TARBALL}.part" "${DB_URL}"
mv "${TARBALL}.part" "${TARBALL}"

echo "[$(date)] extracting"
tar -xzf "${TARBALL}" -C "${KRAKEN_DB}"
rm -f "${TARBALL}"

for f in hash.k2d opts.k2d taxo.k2d; do
    [[ -s "${KRAKEN_DB}/${f}" ]] || die "missing ${f} after extraction — download may be incomplete"
done

echo "[$(date)] DONE. Database in ${KRAKEN_DB}:"
du -sh "${KRAKEN_DB}"
ls -la "${KRAKEN_DB}"
