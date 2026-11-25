library(tidyverse)
library(ggtree)
library(ape)
library(ggplot2)

# ── Paths ─────────────────────────────────────────────────────────────────────
tree_file <- "06_iqtree/gilli_phylo_v2.treefile"
md_host   <- read.table("~/Desktop/Ch2/American_bumblebe_md_2025.tsv",
                        sep = "\t", header = TRUE) %>%
  mutate(across(c(SampleID, SpeciesID, Country), trimws))
bee_pal_raw <- read.csv("~/Desktop/Ch2/paleta_abejorros.csv")
bee_pal     <- setNames(bee_pal_raw$color_hex, bee_pal_raw$especie)
ref_meta    <- read.delim("00_metadata/genome_metadata.tsv", stringsAsFactors = FALSE)
# ANI cluster / GTDB species per MAG, from the per-tip annotation table
# (09_phylogeny/00_metadata/mag_tree_annotation.tsv, rows tree == "gilli_v2")
annot <- read.delim("../00_metadata/mag_tree_annotation.tsv", stringsAsFactors = FALSE) %>%
  filter(tree == "gilli_v2", ani_cluster != "", !ani_cluster %in% c("reference", "outgroup")) %>%
  distinct(mag_id, ani_cluster) %>%
  mutate(ani_cluster = sub("Gilliamella sp. (sin especie GTDB)", "Gilliamella sp. (no GTDB species)", ani_cluster, fixed = TRUE))

# ── Species abbreviations ─────────────────────────────────────────────────────
sp_abbr <- c(
  atratus = "atra", bifarius = "bifa", bifarius_nearcticus = "bifa.near",
  fervidus = "ferv", flavifrons_dimidiatus = "flav.dimi",
  flavifrons_flavifrons = "flav.flav", funebris = "fune",
  griseocollis = "gris", hortulanus = "hort", huntii = "hunt",
  insularis = "insu", mixtus = "mixt", morrisoni = "morr",
  nevadensis = "neva", pensylvanicus = "pens", robustus = "robu",
  rubicundus = "rubi", rufocinctus = "rufo", sonorus = "sono",
  sylvicola = "sylv", vandykei = "vand", vosnesenskii = "vosn",
  Unknown = "unk"
)
country_abbr   <- c(Colombia = "COL", US = "US", Unknown = "unk")
country_shapes <- c(Colombia = 21, US = 24, Unknown = 23, Reference = 22, Outgroup = 25)

# ── Colour palette ────────────────────────────────────────────────────────────
# Gilliamella: refs are mostly Bombus; two strains from stingless bees (ESL0721/0727)
extra_pal <- c(
  "Apis"         = "#E8A823",
  "BombusRef"    = "#9966CC",   # purple for Bombus reference isolates
  "StinglessBee" = "#E84646",   # red for stingless bee Orbaceae (ESL0721/0727)
  "Outgroup"     = "#FFFFFF"   # white: nevadensis already uses grey
)
full_pal <- c(bee_pal, extra_pal)

# ── Ref + outgroup tip metadata ───────────────────────────────────────────────
ref_tips <- ref_meta %>%
  filter(type %in% c("REF", "OUTGROUP")) %>%
  transmute(
    label     = genome_id,   # OrthoFinder uses filename (genome_id) as tip label
    SpeciesID = case_when(
      type == "OUTGROUP"              ~ "Outgroup",
      host_group == "Apis"           ~ "Apis",
      host_group == "Stingless_bee"  ~ "StinglessBee",
      TRUE                           ~ "BombusRef"
    ),
    Country   = ifelse(type == "OUTGROUP", "Outgroup", "Reference"),
    tip_label = paste0(
      genome_id, "\n[",
      ifelse(host_bee == "Unknown" | is.na(host_bee), gill_species, host_bee),
      "]"
    ),
    is_ref    = TRUE
  )

# ── Load + root tree ──────────────────────────────────────────────────────────
tree    <- read.tree(tree_file)

# Root at Pasteurella dagmatis, the only true outgroup (Pasteurellaceae).
# Gilli_Pra-s* are Gilliamella references from Bombus pratorum (China); they were
# typed OUTGROUP in the first version of genome_metadata.tsv and are now REF.
pasteurella_tip <- "Pasteurella_dagmatis"
tree_rooted     <- root(tree, outgroup = pasteurella_tip, resolve.root = TRUE)

tree_sc <- tree_rooted
tree_sc$edge.length <- tree_sc$edge.length ^ 0.35

# UFboot ≥ 70
tree_sc$node.label <- sapply(tree_sc$node.label, function(x) {
  if (is.na(x) || x == "") return("")
  uf <- suppressWarnings(as.numeric(sub(".*/", "", x)))
  if (!is.na(uf) && uf >= 70) as.character(round(uf)) else ""
})

# ── Tip metadata for MAGs ─────────────────────────────────────────────────────
mag_tips <- tibble(label = tree_rooted$tip.label) %>%
  filter(!label %in% ref_tips$label) %>%
  mutate(sample_id = sub("__.*", "", label)) %>%
  left_join(md_host %>% select(SampleID, SpeciesID, Country),
            by = c("sample_id" = "SampleID")) %>%
  mutate(
    SpeciesID = ifelse(is.na(SpeciesID) | SpeciesID == "", "Unknown", SpeciesID),
    Country   = ifelse(is.na(Country)   | Country   == "", "Unknown", Country),
    sp_short  = sp_abbr[SpeciesID],
    ct_short  = country_abbr[Country],
    is_ref    = FALSE
  ) %>%
  group_by(sp_short, ct_short) %>%
  mutate(n = row_number(),
         tip_label = paste0(sp_short, "_", ct_short, "_", sprintf("%02d", n))) %>%
  ungroup() %>%
  left_join(annot, by = c("label" = "mag_id")) %>%
  mutate(tip_label = paste0(tip_label, "  [", ifelse(is.na(ani_cluster), "?", ani_cluster), "]")) %>%
  select(label, SpeciesID, Country, tip_label, is_ref)

all_tips <- bind_rows(mag_tips,
                      ref_tips %>% select(label, SpeciesID, Country, tip_label, is_ref))

# ── Plot: rectangular ─────────────────────────────────────────────────────────
n_taxa <- length(tree_rooted$tip.label)
N_CORE_OGS <- 484   # OrthoFinder core OGs in the v2 concatenation (149 136 aa; verified on the cluster 2026-11-25)
n_mags <- sum(!all_tips$is_ref); n_out <- sum(all_tips$Country == "Outgroup"); n_refs <- sum(all_tips$Country == "Reference")
subtitle_txt <- sprintf("%d core genes · %d taxa (%d MAGs + %d refs + %d outgroup) · IQ-TREE2 UFboot · branches ^0.35",
                        N_CORE_OGS, n_taxa, n_mags, n_refs, n_out)

p_rect <- ggtree(tree_sc, linewidth = 0.35) %<+% all_tips +

  geom_tippoint(aes(fill = SpeciesID, shape = Country),
                size = 3, stroke = 0.4, color = "grey20") +

  geom_tiplab(aes(label = tip_label), size = 2.6, offset = 0.008,
              hjust = 0, repel = TRUE, max.overlaps = Inf) +

  geom_nodelab(size = 2, hjust = 1.2, color = "grey40") +

  geom_treescale(x = 0, y = -1, fontsize = 3, linesize = 0.5,
                 label = "aa subs/site") +

  scale_fill_manual(values = full_pal, name = "Host species / group",
                    na.value = "grey70") +
  scale_shape_manual(values = country_shapes,
                     name   = "Country / type",
                     labels = c(Colombia  = "Colombia (MAG)",
                                US        = "USA (MAG)",
                                Unknown   = "Unknown (MAG)",
                                Reference = "Reference genome",
                                Outgroup  = "Outgroup")) +

  guides(fill  = guide_legend(override.aes = list(shape = 21), ncol = 1),
         shape = guide_legend(override.aes = list(fill  = "grey50"))) +

  labs(title    = expression(italic("Gilliamella") ~ "phylogeny v2 (rooted)"),
       subtitle = subtitle_txt) +

  theme_tree2(legend.position = "right") +
  xlim(0, max(node.depth.edgelength(tree_sc)) * 1.65)

# ── Plot: circular ────────────────────────────────────────────────────────────
p_circ <- ggtree(tree_sc, layout = "circular", linewidth = 0.35) %<+% all_tips +

  geom_tippoint(aes(fill = SpeciesID, shape = Country),
                size = 3, stroke = 0.4, color = "grey20") +

  geom_tiplab(aes(label = tip_label), size = 2.2, offset = 0.018, align = FALSE) +

  scale_fill_manual(values = full_pal, name = "Host species / group",
                    na.value = "grey70") +
  scale_shape_manual(values = country_shapes,
                     name   = "Country / type",
                     labels = c(Colombia  = "Colombia (MAG)",
                                US        = "USA (MAG)",
                                Unknown   = "Unknown (MAG)",
                                Reference = "Reference genome",
                                Outgroup  = "Outgroup")) +

  guides(fill  = guide_legend(override.aes = list(shape = 21), ncol = 1),
         shape = guide_legend(override.aes = list(fill  = "grey50"))) +

  labs(title    = expression(italic("Gilliamella") ~ "phylogeny v2 (circular, rooted)"),
       subtitle = subtitle_txt) +

  theme_tree2(legend.position = "right")

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave("06_iqtree/gilli_phylo_v2_rect.pdf",  p_rect, width = 13, height = 15, device = "pdf")
ggsave("06_iqtree/gilli_phylo_v2_circ.pdf",  p_circ, width = 13, height = 13, device = "pdf")
ggsave("06_iqtree/gilli_phylo_v2_rect.png",  p_rect, width = 13, height = 15, dpi = 200)
ggsave("06_iqtree/gilli_phylo_v2_circ.png",  p_circ, width = 13, height = 13, dpi = 200)

message("Saved: gilli_phylo_v2_rect.pdf/png  &  gilli_phylo_v2_circ.pdf/png")
