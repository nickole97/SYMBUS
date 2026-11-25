#!/bin/bash
# ============================================================
# config.sh
# Author: Nickole Villabona
# Shared configuration for the multi-sample binning experiment
# ============================================================
# Single place to edit paths, the subset, and run options. Every
# script in this directory sources this file.
#
# Usage inside a script:
#   source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
#
# NOT an sbatch script — do not submit it.
# ============================================================

# ======== EXPERIMENT DESIGN ========
# Three arms are compared on the same 22 samples:
#   single-BROKEN : existing results, reused as-is. Reference only, NOT the
#                   baseline — its MaxBin2 was fed the jgi depth file, which
#                   MaxBin2 reads as (contig, contigLen) instead of
#                   (contig, abundance). See ../08_annotation/../../README.
#   single-FIXED  : BASELINE. Single-sample depth, MaxBin2 fed a correct
#                   2-column abundance file. MetaBAT2 and CONCOCT bins are
#                   reused from the existing run — they received correct
#                   input originally, so only MaxBin2 + DAS Tool are rerun.
#   multi         : depth matrix built from all 22 samples per assembly.
#
# The single-sample side is FROZEN. Nothing here writes into
# 03_assembly/, 04_mapping/, 05_binning/, 06_MAGs_selected/ or 07_dRep_*/.

# ======== PATHS: existing pipeline (READ-ONLY) ========
WORK_BASE=/dfs10/hammert-lab/nvillabo/Ch2/work_meta
READS_DIR=${WORK_BASE}/preprocessing/02_host_removal
ASSEMBLY_DIR=${WORK_BASE}/03_assembly
SINGLE_MAP_DIR=${WORK_BASE}/04_mapping
SINGLE_BIN_DIR=${WORK_BASE}/05_binning

# ======== PATHS: this experiment (WRITE) ========
# Bulky intermediates (cross-mapped BAMs, minimap2 indices). ON dfs10, NOT /pub:
# /pub/nvillabo is not this project's working space and holds other projects
# under a shared per-user quota (it overflowed on 2026-11-25 and the user asked
# that it be left alone). 02_depth_matrix deletes each assembly's BAMs as soon
# as its depth products are written, so this directory stays small.
SCRATCH=${WORK_BASE}/12_multisample_test/scratch
# Bins, quality tables and summaries stay on the lab filesystem.
RESULTS=${WORK_BASE}/12_multisample_test

# Depth products (jgi depth matrix, per-sample abundance files, CONCOCT cut-up
# FASTA + coverage table) live on the lab filesystem, NOT on scratch: on
# 2026-11-25 the /pub per-user quota (~75 GB) filled with the 484 BAMs and 02
# could not even write a text file there. Only BAMs and minimap2 indices stay
# on scratch, and 02 deletes each assembly's BAMs as soon as it is done.
DEPTH_DIR=${RESULTS}/02_depth

PAIRS_FILE=${RESULTS}/backmap_pairs.tsv
LOG_DIR=${RESULTS}/logs

# ======== SUBSET ========
# SCRIPT_DIR is set by every sbatch script before sourcing this file (SLURM runs
# scripts from a spool dir, so BASH_SOURCE is unreliable there); fall back to
# BASH_SOURCE only for interactive use on the login node.
SUBSET_TSV="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}/subset_22.tsv"

# Sample IDs only, in a stable order. SLURM array indices are 1-based and
# index THIS list, so never reorder it mid-experiment.
subset_samples() {
    awk -F'\t' 'NR>1 && $1!="" {print $1}' "${SUBSET_TSV}"
}
# Tolerant of a missing TSV so that check_subset can report it properly
# instead of the caller dying inside this assignment under `set -e`.
N_SUBSET=$( { subset_samples 2>/dev/null || true; } | wc -l | tr -d ' ')

# ======== SAMPLE SELECTION (pluggable) ========
# Which selector produces backmap_pairs.tsv. The selector is the ONLY thing
# that changes when scaling to all 113 samples:
#   all_vs_all : every subset sample's reads vs every subset assembly (N x N)
#   simka      : the K most similar samples per assembly (not implemented yet)
SELECTOR=all_vs_all
SIMKA_TOP_K=50          # only used by the simka selector

# ======== RUN OPTIONS ========
ACCOUNT=hammert_lab     # MANDATORY — the personal nvillabo account bills lower
PARTITION=standard
THROTTLE=20             # max concurrent array tasks (sbatch --array=1-N%THROTTLE)

# Delete each cross-mapped BAM as soon as the depth matrix step has consumed
# and validated it. Keeps peak scratch at ~1 assembly's worth of BAMs (~3 GB)
# instead of all 484 (~70 GB). Set to 0 to keep every BAM.
DELETE_BAMS_AFTER_DEPTH=1

# ======== PYTHON HYGIENE ========
# The user's ~/.local/lib/python3.10/site-packages holds numpy 2.2.6 (pip
# --user), which the env's python loads BEFORE the env's own numpy 1.26.4;
# CONCOCT's compiled vbgmm then dies with "compiled using NumPy 1.x cannot be
# run in NumPy 2.2.6" (03c, 2026-11-25). Disable user site-packages for every
# step-12 job. Do not "fix" this by touching ~/.local.
export PYTHONNOUSERSITE=1

# ======== TOOLS ========
# Same absolute paths the existing pipeline uses, so this experiment runs the
# identical binner builds. Changing a binner version would confound the test.
CONDA_BIN=/data/homezvol3/nvillabo/.conda/envs/binning_env/bin
METABAT2=${CONDA_BIN}/metabat2
MAXBIN2=${CONDA_BIN}/run_MaxBin.pl

# jgi_summarize_bam_contig_depths exists in BOTH binning_env (ships MetaBAT
# 2.18, 2025) and quast_env (ships MetaBAT 2.12.1, 2017), and the two binaries
# differ. Step 3a of the original pipeline generated every existing depth file
# with the quast_env copy, so we use that one here too: the single-FIXED arm
# reuses those old depth files, and switching generators for the multi arm
# would add a second variable to the comparison. Do not "upgrade" this path.
JGI_DEPTH=/data/homezvol3/nvillabo/.conda/envs/quast_env/bin/jgi_summarize_bam_contig_depths
CONCOCT_CUTUP=${CONDA_BIN}/cut_up_fasta.py
CONCOCT_COVTAB=${CONDA_BIN}/concoct_coverage_table.py
CONCOCT=${CONDA_BIN}/concoct
CONCOCT_MERGE=${CONDA_BIN}/merge_cutup_clustering.py
CONCOCT_EXTRACT=${CONDA_BIN}/extract_fasta_bins.py
FASTA2BIN=${CONDA_BIN}/Fasta_to_Contig2Bin.sh
DASTOOL=${CONDA_BIN}/DAS_Tool

CHECKM=/data/homezvol3/nvillabo/.conda/envs/checkm_env/bin/checkm
CHECKM_DATA=/data/homezvol3/nvillabo/.conda/envs/checkm_env/checkm_data
# CheckM1 on every arm, deliberately. CheckM2 is the better tool but mixing
# versions across arms would invalidate the comparison.

GTDBTK_ENV=gtdbtk-2.5.2

# ======== BINNER PARAMETERS (must match the single-sample run) ========
MIN_CONTIG_LEN=1500     # MetaBAT2 -m, MaxBin2 -min_contig_length
CONCOCT_CHUNK=10000     # cut_up_fasta.py -c

# ======== HELPERS ========
# Resolve a sample's array index (1-based) or its ID from an index.
sample_from_index() { subset_samples | sed -n "${1}p"; }

die() { echo "[ERROR] $*" >&2; exit 1; }

check_subset() {
    [[ -s "${SUBSET_TSV}" ]] || die "subset TSV not found: ${SUBSET_TSV}"
    [[ "${N_SUBSET}" -gt 0 ]] || die "subset is empty: ${SUBSET_TSV}"
}
