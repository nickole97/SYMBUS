library(tidyverse)
library(ggtree)
library(ape)
library(ggplot2)

# ── Paths ─────────────────────────────────────────────────────────────────────
tree_file <- "03_iqtree/snod_phylo_v4.treefile"
md_host   <- read.table("~/Desktop/Ch2/American_bumblebe_md_2025.tsv",
                        sep = "\t", header = TRUE) %>%
  mutate(across(c(SampleID, SpeciesID, Country), trimws))
bee_pal_raw <- read.csv("~/Desktop/Ch2/paleta_abejorros.csv")
bee_pal     <- setNames(bee_pal_raw$color_hex, bee_pal_raw$especie)

# v4 reuses v2 ingroup metadata + v3 outgroup rows
ref_meta_v2 <- read.delim("../snod_phylo_v2/00_metadata/genome_metadata.tsv",
                           stringsAsFactors = FALSE)
ref_meta_v3 <- read.delim("../snod_phylo_v3/00_metadata/genome_metadata.tsv",
                           stringsAsFactors = FALSE) %>%
  filter(type == "OUTGROUP")
ref_meta <- bind_rows(ref_meta_v2, ref_meta_v3)

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

# ── Extended colour palette ───────────────────────────────────────────────────
extra_pal <- c(
  "Apis"         = "#E8A823",
  "StinglessBee" = "#E84646",
  "pascuorum"    = "#7FBFFF",
  "Outgroup"     = "#BBBBBB"
)
full_pal <- c(bee_pal, extra_pal)

# ── Ref + outgroup tip metadata ───────────────────────────────────────────────
ref_tips <- ref_meta %>%
  filter(type %in% c("REF", "OUTGROUP")) %>%
  transmute(
    label     = accession,
    SpeciesID = case_when(
      type == "OUTGROUP"              ~ "Outgroup",
      host_group == "Apis"           ~ "Apis",
      host_group == "Stingless_bee"  ~ "StinglessBee",
      genome_id  == "S_gandavensis_type" ~ "pascuorum",
      genome_id  == "S_communis_type"    ~ "griseocollis",
      TRUE ~ "Unknown"
    ),
    Country   = ifelse(type == "OUTGROUP", "Outgroup", "Reference"),
    tip_label = paste0(snod_species, "\n[", host_bee, "]"),
    is_ref    = TRUE
  )

# ── Load + root tree ──────────────────────────────────────────────────────────
tree <- read.tree(tree_file)

# Root at LCA of outgroup taxa (Simonsiella + Rivihabitans)
outgroup_tips <- ref_meta %>% filter(type == "OUTGROUP") %>% pull(accession)
outgroup_node <- getMRCA(tree, outgroup_tips)
tree_rooted   <- root(tree, node = outgroup_node, resolve.root = TRUE)

tree_sc <- tree_rooted
tree_sc$edge.length <- tree_sc$edge.length ^ 0.35

# UFboot ≥ 70 on internal nodes
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
  select(label, SpeciesID, Country, tip_label, is_ref)

n_mags <- nrow(mag_tips)
n_refs <- ref_tips %>% filter(is_ref, !grepl("Outgroup", Country)) %>% nrow()
n_out  <- ref_tips %>% filter(!is_ref | grepl("Outgroup", Country)) %>% nrow()
all_tips <- bind_rows(mag_tips,
                      ref_tips %>% select(label, SpeciesID, Country, tip_label, is_ref))

# ── Plot: rectangular ─────────────────────────────────────────────────────────
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

  labs(
    title    = expression(italic("Snodgrassella") ~ "phylogeny v4 (rooted, best-practice)"),
    subtitle = paste0(
      "411 core genes · ", nrow(all_tips), " taxa (",
      n_mags, " MAGs + ", n_refs, " refs + 2 outgroups) · ",
      "ingroup-only OrthoFinder + BLAST outgroups · IQ-TREE2 UFboot · branches ^0.35"
    )
  ) +

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

  labs(title    = expression(italic("Snodgrassella") ~ "phylogeny v4 (circular, rooted)"),
       subtitle = "411 core genes · ingroup-only OrthoFinder + BLAST outgroups · IQ-TREE2 UFboot · branches ^0.35") +

  theme_tree2(legend.position = "right")

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave("03_iqtree/snod_phylo_v4_rect.pdf",  p_rect, width = 13, height = 15, device = "pdf")
ggsave("03_iqtree/snod_phylo_v4_circ.pdf",  p_circ, width = 13, height = 13, device = "pdf")
ggsave("03_iqtree/snod_phylo_v4_rect.png",  p_rect, width = 13, height = 15, dpi = 200)
ggsave("03_iqtree/snod_phylo_v4_circ.png",  p_circ, width = 13, height = 13, dpi = 200)

message("Saved: snod_phylo_v4_rect.pdf/png  &  snod_phylo_v4_circ.pdf/png")
