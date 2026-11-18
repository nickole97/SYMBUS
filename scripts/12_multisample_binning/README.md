# Step 12 — Multi-sample binning experiment

A **test**, not a pipeline change. Steps 01–11 are untouched and the
single-sample results are frozen: nothing here writes into `03_assembly/`,
`04_mapping/`, `05_binning/`, `06_MAGs_selected/` or `07_dRep_*/`.

Suggested by Charlotte Huyghe (Engel lab, Lausanne): keep the per-sample
MEGAHIT assemblies exactly as they are, but build the binning depth matrix by
backmapping **many** samples' reads to each assembly, so coverage covariance
across samples gives the binners more signal. This is **not** co-assembly.

Supporting benchmark: Han H, Wang Z, Zhu S (2025) *Benchmarking metagenomic
binning tools on real datasets across sequencing platforms and binning modes.*
Nat Commun 16:2865 — benchmarks this exact binner stack (CONCOCT, MaxBin2,
MetaBAT2, DAS Tool) and reports +100% medium-quality and +194% near-complete
MAGs for multi-sample over single-sample on short-read data, with gains
becoming significant at 15+ samples.

## Why there are three arms

| Arm | Depth | MaxBin2 abundance | Runs needed |
|---|---|---|---|
| single-BROKEN | 1 sample | jgi table → contigLen read as abundance | none; reused as a reference point |
| **single-FIXED** ← baseline | 1 sample | correct 2-column file | MaxBin2 + DAS Tool only |
| multi | 22 samples | correct per-sample files via `-abund_list` | cross-map + 3 binners + DAS Tool |

MaxBin2 supplies **18.6%** of the DAS Tool bins dataset-wide (34 of 183, 16 of
them high quality, and the *only* source of bins in 8 samples). Fixing its
abundance format in the multi arm while leaving it broken in the single arm
would let the bug fix masquerade as a method effect — hence the baseline is
single-FIXED, and single-BROKEN is kept only for reference.

**single-FIXED is cheap** because MetaBAT2 and CONCOCT received correct input in
the original run. Only MaxBin2 is rebuilt; `04b` reuses the existing MetaBAT2
and CONCOCT bins straight from `05_binning/`, read-only.

## Subset: 22 samples

Defined in `subset_22.tsv`, which also freezes each sample's pre-experiment bin
count so the report can compare against it. **Array indices index that file, so
never reorder it mid-experiment.**

| Role | n | What |
|---|---|---|
| `failure` | 14 | *Pyrobombus* with 0–3 bins and 0.6–2.8M non-host reads (flavifrons ×4, huntii ×4, bifarius ×3, vandykei ×2, mixtus ×1) |
| `control` | 6 | *vosnesenskii*, 2–5 bins each — tests whether multi-sample depth **degrades** samples that already work |
| `depth_contrast` | 2 | *rufocinctus* (Cullumanobombus), 4.7–5.7M reads yet zero bins — a failure mode that is not depth |

Only 7 *Pyrobombus* samples dataset-wide have ≥2 bins and 6 are *vosnesenskii*,
so within-subgenus positive controls could not be anything else. Baseline across
the subset: 25 bins.

## Order of operations

```bash
# 0. pair list (login node, instant)
bash 00_make_backmap_pairs.sh

# 1-2. cross-map, then derive every coverage format and release the BAMs
sbatch --array=1-484%20 01_crossmap.sh
sbatch --array=1-22%20  02_depth_matrix.sh

# 3. binning — multi arm (a,b,c) and the fixed single-sample arm (d)
sbatch --array=1-22%20 03a_metabat_multi.sh
sbatch --array=1-22%20 03b_maxbin_multi.sh
sbatch --array=1-22%20 03c_concoct_multi.sh
sbatch --array=1-22%20 03d_maxbin_single_fixed.sh

# 4. refinement, both arms
sbatch --array=1-22%20 04a_dastool_multi.sh
sbatch --array=1-22%20 04b_dastool_single_fixed.sh

# 5-6. quality and taxonomy
sbatch --array=1-2 05_checkm1_both.sh
sbatch           06_gtdbtk_both.sh

# 7. report
Rscript -e 'rmarkdown::render("07_compare_report.Rmd")'

# 8. independent: does the rufocinctus contrast mean what it looks like?
sbatch 08a_kraken2_setup.sh      # once, ~8 GB download
sbatch --array=1-4 08b_kraken2_screen.sh
```

`00` refuses to regenerate the pair list once `01` has started: array indices
point at lines in that file, and regenerating it would silently remap task
numbers to different pairs.

## Pluggable sample selection

The seam is a file, not a function. The selector writes
`backmap_pairs.tsv` with two columns — `assembly_sample`, `read_sample` — and
nothing downstream cares how it was produced.

- `selectors/all_vs_all.sh` — 22 × 22 = 484 pairs. Used now; Simka is pointless
  at this scale.
- `selectors/simka.sh` — **stub**. For all 113 samples, where all-vs-all would be
  12,769 mappings. Charlotte's pipeline backmaps the ~50 most similar samples
  per assembly, selected with Simka. The header documents the interface.

Switching is one line: `SELECTOR=` in `config.sh`.

## Where things go

| | Path |
|---|---|
| Bulky intermediates (BAMs, depth matrices, cut-up FASTAs) | `/pub/nvillabo/Ch2/12_multisample_test/` (scratch) |
| Bins, CheckM, GTDB-Tk, summaries | `work_meta/12_multisample_test/` |
| Scripts and the report source | this directory |

`02` deletes each assembly's BAMs as soon as all three coverage products are
built and validated (`DELETE_BAMS_AFTER_DEPTH=1`), keeping peak scratch at about
one assembly's worth (~3 GB) instead of all 484 (~70 GB). Set it to `0` to keep
them.

## Traps worth knowing about

- **`02` is the only step that reads the BAMs.** The CONCOCT coverage table is
  built there rather than in `03c` precisely so the BAMs can be deleted
  afterwards. Do not move it.
- **The jgi binary is pinned to `quast_env`.** `jgi_summarize_bam_contig_depths`
  exists in both `binning_env` (MetaBAT 2.18, 2025) and `quast_env`
  (MetaBAT 2.12.1, 2017), and they differ. Every existing depth file was made
  with the `quast_env` copy, and the single-FIXED arm reuses those files, so the
  multi arm must use the same binary or the comparison gains a second variable.
  `config.sh` says *do not upgrade this path* for this reason.
- **`03b` refuses to run** unless every abundance file has exactly 2 columns.
  The original bug was completely silent: MaxBin2 does not complain, it just
  uses contig length as abundance.
- **Conda before `set -euo pipefail`** in `04a`, `04b` and `06` — activation
  fails under SLURM otherwise.
- **CheckM1, not CheckM2**, on every arm. CheckM2 is the better tool but mixing
  versions across arms would invalidate the comparison.
- **256 GB for `06`** is deliberate. A validation run peaked at ~110 GiB, but
  both its genomes were resolved by the skani ANI pre-screen, which skips
  pplacer. Novel bins — what this experiment hopes to recover — take the pplacer
  path instead.

## What the report answers

1. Bins per sample, three arms, paired.
2. MIMAG tiers (CheckM1): HQ >90%/<5%, MQ ≥50%/<10%.
3. **Per-binner contribution to DAS Tool output, single-FIXED vs multi.** CONCOCT
   is *already* the largest contributor (86 of 183 bins, 59 high quality); the
   hypothesis is that its **retention rate** improves — it currently produces
   16–55 bins per sample of which DAS Tool keeps 0–2.
4. **Taxa recovered** — set difference of GTDB labels at genus and species level,
   so extra copies of an already-recovered species do not count as a gain. This
   is the metric the experiment exists for.
5. **Degradation** of the 6 controls, which decides whether rerunning all 113 is
   a good trade. Caveat stated in the report: the controls are from a different
   collection than the failure block, so within the all-vs-all they covary mostly
   with each other and receive a weaker perturbation.
6. Within-subgenus vs cross-subgenus (20 *Pyrobombus* vs 2 *Cullumanobombus*).
   Charlotte's view is that cross-sample noise is not a concern and is partly
   necessary; this tests that on our own data.
