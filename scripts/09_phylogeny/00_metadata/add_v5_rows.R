#!/usr/bin/env Rscript
# add_v5_rows.R — append the snod_v5 rows (62 tips) to mag_tree_annotation.tsv.
# MAG rows copy the sample/host/GTDB/CheckM fields from the snod_v3 rows (same
# MAG) or from the "none" rows (the 6 MAGs that were in no tree before v5);
# the tree-specific fields (tree, tree_file, label_in_plot, notes) come from
# the v5 tip table. Reference and outgroup rows come from v5 genome_metadata.
# Idempotent: existing snod_v5 rows are dropped before appending.
suppressMessages(library(dplyr))
here <- "/Users/nickolevillabona/Desktop/Ch2/SYMBUS/scripts/09_phylogeny"
ann  <- read.delim(file.path(here, "00_metadata/mag_tree_annotation.tsv"), stringsAsFactors = FALSE, na.strings = "", colClasses = "character")
tips <- read.delim(file.path(here, "snod_phylo_v5/07_iqtree/snod_phylo_v5_tip_table.tsv"), stringsAsFactors = FALSE)
gm   <- read.delim(file.path(here, "snod_phylo_v5/00_metadata/genome_metadata.tsv"), stringsAsFactors = FALSE)
ann  <- ann %>% filter(tree != "snod_v5")
tree_file <- "SYMBUS/scripts/09_phylogeny/snod_phylo_v5/07_iqtree/snod_phylo_v5.treefile (cluster 2026-11-25; 533 OGs, 127575 aa, 62 taxa)"
# source row per MAG: prefer snod_v3, else "none"
src <- ann %>% filter(tree %in% c("snod_v3", "none")) %>%
  arrange(mag_id, factor(tree, levels = c("snod_v3", "none"))) %>% distinct(mag_id, .keep_all = TRUE)
mag_tips <- tips %>% filter(!grepl("^(GC[AF]_|Simonsiella|Rivihabitans)", label))
stopifnot(all(mag_tips$label %in% src$mag_id))
mag_rows <- src %>% filter(mag_id %in% mag_tips$label) %>%
  mutate(tree = "snod_v5", tree_file = tree_file, tip_label_in_tree = mag_id,
         label_in_plot = tips$tip_code[match(mag_id, tips$label)],
         ani_cluster   = tips$ani_cluster[match(mag_id, tips$label)],
         in_prelim_bac120 = ifelse(is.na(in_prelim_bac120), "", in_prelim_bac120),
         notes = paste0("label in plot = ", label_in_plot,
                        ifelse(tree == "snod_v5" & mag_id %in% ann$mag_id[ann$tree == "none"],
                               "; first tree that includes this MAG (absent from v2/v3)", "")))
ref_tips <- tips %>% filter(grepl("^(GC[AF]_|Simonsiella|Rivihabitans)", label))
# REF tips are labelled by accession, OUTGROUP tips by genome_id: key gm by both
gm <- bind_rows(gm %>% mutate(genome_id = ifelse(type == "REF", accession, genome_id)), gm) %>% distinct(genome_id, .keep_all = TRUE)
ref_rows <- ref_tips %>% transmute(
  tip_label_in_tree = label, tree = "snod_v5", tree_file = tree_file, mag_id = label, sample_id = NA,
  host_species = gm$host_bee[match(label, gm$genome_id)], host_subgenus = NA,
  country = gm$country[match(label, gm$genome_id)],
  gtdb_classification = gm$snod_species[match(label, gm$genome_id)], fastani_ref_ani = NA,
  checkm_completeness = NA, checkm_contamination = NA,
  ani_cluster = ifelse(gm$type[match(label, gm$genome_id)] == "OUTGROUP", "outgroup", "reference"),
  reference_accession = gm$accession[match(label, gm$genome_id)], in_prelim_bac120 = "",
  label_in_plot = label, notes = paste0(gm$type[match(label, gm$genome_id)], "; ", gm$source[match(label, gm$genome_id)], "; ", gm$notes[match(label, gm$genome_id)]))
stopifnot(!any(is.na(ref_rows$host_species)))
out <- bind_rows(ann, mag_rows[names(ann)], ref_rows[names(ann)])
write.table(out, file.path(here, "00_metadata/mag_tree_annotation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat("snod_v5 rows added:", nrow(mag_rows) + nrow(ref_rows), "(MAGs", nrow(mag_rows), ", refs/outgroups", nrow(ref_rows), "); total rows:", nrow(out), "\n")
print(table(out$tree))
