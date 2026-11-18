#!/bin/bash
# ============================================================
# selectors/simka.sh
# Author: Nickole Villabona
# Sample selector: the K most similar samples per assembly (NOT IMPLEMENTED)
# ============================================================
# Placeholder documenting the interface, for scaling this experiment from
# the 22-sample subset to all 113 samples.
#
# Why it is needed at 113 and not at 22:
#   all_vs_all at 113 samples = 12,769 mappings. Infeasible, and mostly
#   wasted: samples sharing no strains with an assembly contribute noise
#   to the coverage covariance rather than signal. Charlotte Huyghe's
#   pipeline (Engel lab) backmaps the ~50 most similar samples per
#   assembly, selected with Simka.
#
# What this script must do when implemented:
#   1. Run Simka once over all samples' non-host reads to get a pairwise
#      k-mer distance matrix (Bray-Curtis on k-mer abundance).
#        simka -in <manifest.txt> -out <dir> -out-tmp <scratch>
#      The manifest is one line per sample: "<id>: <R1>;<R2>"
#   2. For each assembly sample, take the SIMKA_TOP_K nearest samples by
#      that distance, always including the assembly's own reads.
#   3. Emit the pairs on stdout in the same format all_vs_all.sh uses:
#        <assembly_sample>\t<read_sample>
#
# Nothing downstream changes: 01_crossmap.sh and 02_depth_matrix.sh read
# backmap_pairs.tsv and do not care how it was produced. Switching is a
# one-line edit of SELECTOR in config.sh.
#
# Open questions to settle before implementing:
#   - Simka on non-host reads or on the trimmed reads? Non-host is the
#     consistent choice, since that is what gets assembled and mapped.
#   - Is K=50 right for this dataset? Charlotte's 50 was tuned on honeybee
#     data with more samples per site. The subset experiment should give a
#     sense of how quickly the gain saturates with the number of samples.
#   - Simka needs its own conda env; it is not in binning_env.
# ============================================================

set -euo pipefail

cat >&2 <<'MSG'
[ERROR] selectors/simka.sh is not implemented yet.

Use SELECTOR=all_vs_all in config.sh for the 22-sample subset.
Implement this selector before scaling to all 113 samples; see the header
comments in this file for the required interface.
MSG
exit 1
