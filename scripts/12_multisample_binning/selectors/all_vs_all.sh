#!/bin/bash
# ============================================================
# selectors/all_vs_all.sh
# Author: Nickole Villabona
# Sample selector: every subset sample against every subset assembly
# ============================================================
# Emits the backmap pair list on stdout: for each assembly, every sample
# in the subset (including itself) contributes its reads.
#
# With 22 samples this is 484 pairs, which is cheap. Simka-based selection
# is unnecessary at this scale — see selectors/simka.sh for the version
# that matters when scaling to all 113 samples.
#
# Output format (stdout, tab-separated, no header):
#   <assembly_sample>\t<read_sample>
#
# Called by 00_make_backmap_pairs.sh. Not submitted directly.
# ============================================================

set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"
check_subset

mapfile -t SAMPLES < <(subset_samples)

for asm in "${SAMPLES[@]}"; do
    for rd in "${SAMPLES[@]}"; do
        printf '%s\t%s\n' "${asm}" "${rd}"
    done
done
