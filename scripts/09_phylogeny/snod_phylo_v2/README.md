# Snodgrassella Phylogeny v2 — Expanded with Apis + Stingless Bee References

## What's new vs. v1 (`snod_phylo/`)

| | v1 | v2 |
|---|---|---|
| Genomes | 7 focal MAGs only | 7 focal MAGs + up to 46 MAGs + references |
| Apis | ❌ | ✅ 5 genomes (S. alvi wkB2, K-FP26, cerana, florea, andreniformis) |
| Named species anchors | ❌ | ✅ S. communis type + S. gandavensis type |
| Stingless bees | ❌ | ✅ ESL0693 (Snodgrassella), ESL0721/ESL0727 (Gilliamella) |
| Host tree | — | YY's mitochondrial tree (TreeFile.nex.con_BombusPhylogeny07) |

## Genomes included

See `00_metadata/genome_metadata.tsv` for the full table.

**Our MAGs (7 focal):**
- Novel sp. A–E (Colombian Bombus, novel Snodgrassella)
- S. alvi_D (202006191520__22_sub)
- S. communis (202206071041G__39)

**Reference genomes:**
- S. alvi wkB2 — type strain (GCF_000600005.1)
- S. alvi K-FP26 — A. mellifera iberiensis (CP132376)
- novel Apis group — A. cerana, A. florea, A. andreniformis
- S. communis type strain — B. griseocollis (GCF_002777335.1)
- S. gandavensis type strain — B. pascuorum (GCA_914768025.1)
- ESL0693 — Snodgrassella from stingless bee (TaxID 2983239)
- ESL0721/ESL0727 — Gilliamella from stingless bee (TaxID 2983240/2983241)

## Pipeline steps

```
01_reference_genomes/   ← downloaded from NCBI (scripts/00_download_refs.sh)
02_our_MAGs/            ← our MAG FASTAs (scripts/00b_link_our_MAGs.sh)
        ↓
03_prodigal/            ← gene prediction (.faa) for all genomes
        ↓
04_orthofinder/         ← single-copy core orthologs
        ↓
05_core_OGs/            ← filtered OGs (all genomes present)
        ↓
06_alignments/          ← MAFFT + trimAl + concatenated supermatrix
        ↓
07_iqtree/              ← ML tree, -m TEST, 1000 bootstrap
```

## Pending before running

- [x] Confirm GCA accessions for ESL0693 (GCA_029433475.1), ESL0721 (GCA_029433365.1), ESL0727 (GCA_029433345.1)
- [x] Scope: all 46 MAGs from original phylogeny (all_snod_mags.txt) + 10 reference genomes = 56 total
- [ ] Verify MAG FASTAs exist at $HOME/Ch2/06_MAGs_selected/ on HPC (check names match all_snod_mags.txt)
- [ ] Install conda environments on HPC (see snod_phylo/install_envs.sh — same envs)
- [ ] Sync this directory to HPC with hpc-repo-sync
