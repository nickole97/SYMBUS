#!/usr/bin/env Rscript
# 07b_quick_compare.R — quick two-arm comparison (single-FIXED vs multi) from the
# CheckM1 tables, the GTDB-Tk summary and the DAS Tool summaries, run locally on
# copies of: ${RESULTS}/checkm/checkm_qa_{single_fixed,multi}.tsv,
# ${RESULTS}/gtdbtk/gtdbtk.bac120.summary.tsv, ${RESULTS}/{single_fixed,multi}/*/DASTool/*_summary.tsv.
# Set S to the folder holding those copies. First run 2026-11-25: same 25 bins and
# taxa in both arms; multi improves quality (HQ 16 vs 15, LQ 2 vs 5; mean
# completeness 87.8 vs 84.5 %). The full report is 07_compare_report.Rmd.
suppressMessages({library(dplyr); library(tidyr); library(stringr)})
S <- "/private/tmp/claude-502/-Users-nickolevillabona-Desktop-Ch2/ba8ff9a8-f13c-4971-beb4-b47090bbdbce/scratchpad/binning_results"
rd <- function(f, arm) read.delim(file.path(S, f), check.names = FALSE) %>% mutate(arm = arm)
ck <- bind_rows(rd("checkm_qa_single_fixed.tsv", "single_fixed"), rd("checkm_qa_multi.tsv", "multi")) %>%
  rename(comp = Completeness, cont = Contamination, size = `Genome size (bp)`, n_scaf = `# scaffolds`) %>%
  mutate(tier = case_when(comp >= 90 & cont < 5 ~ "HQ", comp >= 50 & cont < 10 ~ "MQ", TRUE ~ "LQ"))
cat("=== bins y calidad por brazo\n")
print(ck %>% group_by(arm) %>% summarise(n = n(), HQ = sum(tier=="HQ"), MQ = sum(tier=="MQ"), LQ = sum(tier=="LQ"),
  mean_comp = round(mean(comp),1), mean_cont = round(mean(cont),2), median_size_Mb = round(median(size)/1e6,2), median_scaf = median(n_scaf)) %>% as.data.frame())
g <- read.delim(file.path(S, "gtdbtk.bac120.summary.tsv")) %>% select(user_genome, classification) %>%
  mutate(arm = ifelse(str_starts(user_genome, "ms__"), "multi", "single_fixed"),
         BIN = str_remove(user_genome, "^(ms|sf)__[^_]+__"), SAMPLE = str_match(user_genome, "^(?:ms|sf)__([^_]+)__")[,2],
         genus = str_extract(classification, "(?<=g__)[^;]+"), species = str_extract(classification, "(?<=s__).+"))
cat("\n=== prefijos GTDB:", paste(unique(str_extract(g$user_genome, "^[a-z]+__")), collapse=", "), "| n por brazo:"); print(table(g$arm))
cat("\n=== taxones por brazo (género: n bins)\n"); print(g %>% count(arm, genus) %>% pivot_wider(names_from = arm, values_from = n, values_fill = 0) %>% as.data.frame())
sp <- g %>% mutate(species = ifelse(is.na(species) | species == "", paste0(genus, " (sin especie)"), species)) %>% distinct(arm, species)
cat("\n=== especies solo en un brazo\n"); print(list(solo_multi = setdiff(sp$species[sp$arm=="multi"], sp$species[sp$arm=="single_fixed"]), solo_single_fixed = setdiff(sp$species[sp$arm=="single_fixed"], sp$species[sp$arm=="multi"])))
# por muestra: emparejar bins por especie y comparar calidad
j <- ck %>% left_join(g %>% select(arm, SAMPLE, BIN, genus, species), by = c("arm","SAMPLE","BIN"))
cat("\n=== bins sin clasificación GTDB:", sum(is.na(j$genus)), "\n")
w <- j %>% mutate(sp = ifelse(is.na(species) | species=="", paste0(genus," sp."), species)) %>% group_by(SAMPLE, sp, arm) %>% summarise(comp = max(comp), cont = min(cont), n=n(), .groups="drop") %>%
  pivot_wider(names_from = arm, values_from = c(comp, cont, n))
cat("\n=== por muestra y especie: completitud / contaminación single_fixed vs multi\n")
print(w %>% mutate(d_comp = round(comp_multi - comp_single_fixed,1), d_cont = round(cont_multi - cont_single_fixed,2)) %>% select(SAMPLE, sp, comp_single_fixed, comp_multi, d_comp, cont_single_fixed, cont_multi, d_cont) %>% arrange(SAMPLE, sp) %>% as.data.frame(), row.names = FALSE)
cat("\n=== resumen de deltas (bins emparejados):"); d <- w %>% filter(!is.na(comp_multi), !is.na(comp_single_fixed)); cat(" n pares =", nrow(d), "| mediana d_comp =", median(d$comp_multi - d$comp_single_fixed), "| mediana d_cont =", median(d$cont_multi - d$cont_single_fixed), "| mejor comp en multi:", sum(d$comp_multi > d$comp_single_fixed + 1), "| peor:", sum(d$comp_multi < d$comp_single_fixed - 1), "\n")
# binner de origen
das <- list.files(file.path(S,"das"), pattern = "summary.tsv$", recursive = TRUE, full.names = TRUE)
bs <- bind_rows(lapply(das, function(f) read.delim(f, colClasses = c(bin = "character")) %>% mutate(arm = ifelse(grepl("/multi/", f), "multi", "single_fixed"), SAMPLE = basename(dirname(dirname(f))))))
cat("\n=== binner de origen por brazo\n"); print(bs %>% count(arm, bin_set) %>% pivot_wider(names_from = arm, values_from = n, values_fill = 0) %>% as.data.frame())
