#!/usr/bin/env Rscript

library(tidyverse)
library(decontam)
library(gridExtra)
library(grid)

# =====================
# PARSE ARGUMENTS
# =====================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: contamination_filtering.R <asv_table> <taxonomy> <kraken2_output> <metadata>")
}

counts_file   <- args[1]
taxonomy_file <- args[2]
kraken_file   <- args[3]
metadata_file <- args[4]

# =====================
# PREPARE DISCARD LISTS
# =====================

alien_taxa <- c(
  "Homo sapiens",
  "Thermus aquaticus"
)

CONTAMINATION_THRESHOLD <- 15.0

# ===============
# LOAD INPUT DATA
# ===============

counts_df   <- read_tsv(counts_file, show_col_types = FALSE) %>%
  rename(ASV_ID = `#OTU ID`)
taxonomy_df <- read_tsv(taxonomy_file, show_col_types = FALSE) %>%
  rename(ASV_ID = 1)
kraken_df   <- read_tsv(kraken_file,
                        col_names = c("Status", "ASV_ID", "Taxonomy", "Length", "LCA"),
                        show_col_types = FALSE)
metadata_df <- read_tsv(metadata_file, show_col_types = FALSE) %>%
  mutate(across(everything(), as.character))

# ========================
# BUILD SAMPLE METADATA
# ========================

sample_metadata <- metadata_df %>%
  select(`sample-id`, healthy) %>%
  mutate(sample_type = case_when(
    healthy == "blank" ~ "Blank",
    healthy == "yes"   ~ "Healthy",
    healthy == "no"    ~ "Non-healthy"
  ))

# =======================
# FILTER STEP 1: DECONTAM
# =======================

asv_matrix <- counts_df %>%
  column_to_rownames("ASV_ID") %>%
  as.matrix() %>%
  t()

sample_names <- rownames(asv_matrix)
is_blank     <- sample_names %in%
  sample_metadata$`sample-id`[sample_metadata$healthy == "blank"]

contamination_results <- isContaminant(asv_matrix, neg = is_blank,
                                       method = "prevalence", threshold = 0.1)

kitome_asvs <- rownames(contamination_results[contamination_results$contaminant == TRUE, ])

cat(sprintf("Decontam identified %d contaminant ASVs.\n", length(kitome_asvs)))

# ======================
# FILTER STEP 2: KRAKEN2
# ======================

alien_asvs <- kraken_df %>%
  filter(str_detect(Taxonomy, paste(alien_taxa, collapse = "|"))) %>%
  pull(ASV_ID)

cat(sprintf("Identified %d biological anomaly ASVs.\n", length(alien_asvs)))

# ==========================
# FLAG CONTAMINANTS
# ==========================

apply_flags <- function(df) {
  df %>%
    mutate(Contamination_Flag = case_when(
      ASV_ID %in% kitome_asvs & ASV_ID %in% alien_asvs ~ "Flagged_Both",
      ASV_ID %in% kitome_asvs                           ~ "Flagged_Kitome",
      ASV_ID %in% alien_asvs                            ~ "Flagged_Alien",
      TRUE                                               ~ "Passed"
    ))
}

annotated_counts_df <- apply_flags(counts_df)

# ==========================
# SUMMARY TABLE PER SAMPLE
# ==========================

counts_long <- annotated_counts_df %>%
  pivot_longer(-c(ASV_ID, Contamination_Flag),
               names_to  = "sample-id",
               values_to = "count") %>%
  filter(count > 0) %>%
  left_join(sample_metadata, by = "sample-id")

contamination_summary <- counts_long %>%
  group_by(`sample-id`, sample_type) %>%
  summarise(
    total_asvs       = n(),
    n_passed         = sum(Contamination_Flag == "Passed"),
    n_flagged_kitome = sum(Contamination_Flag == "Flagged_Kitome"),
    n_flagged_alien  = sum(Contamination_Flag == "Flagged_Alien"),
    n_flagged_both   = sum(Contamination_Flag == "Flagged_Both"),
    pct_flagged      = round(100 * (1 - n_passed / total_asvs), 2),
    .groups = "drop"
  )

# ==========================
# TABLE 1 — GROUP AVERAGES
# (Healthy and Non-healthy only)
# ==========================

group_averages <- contamination_summary %>%
  filter(sample_type %in% c("Healthy", "Non-healthy")) %>%
  group_by(sample_type) %>%
  summarise(
    total_asvs       = round(mean(total_asvs), 1),
    n_passed         = round(mean(n_passed), 1),
    n_flagged_kitome = round(mean(n_flagged_kitome), 1),
    n_flagged_alien  = round(mean(n_flagged_alien), 1),
    n_flagged_both   = round(mean(n_flagged_both), 1),
    pct_flagged      = round(mean(pct_flagged), 2),
    .groups = "drop"
  ) %>%
  rename(
    "Sample Type"    = sample_type,
    "Avg Total ASVs" = total_asvs,
    "Avg Passed"     = n_passed,
    "Avg Kitome"     = n_flagged_kitome,
    "Avg Alien"      = n_flagged_alien,
    "Avg Both"       = n_flagged_both,
    "Avg % Flagged"  = pct_flagged
  )

# ==========================
# TABLE 2 — HIGH CONTAMINATION SAMPLES
# (Healthy and Non-healthy with pct_flagged > threshold)
# ==========================

high_contam <- contamination_summary %>%
  filter(
    sample_type %in% c("Healthy", "Non-healthy"),
    pct_flagged > CONTAMINATION_THRESHOLD
  ) %>%
  arrange(desc(pct_flagged)) %>%
  select(`sample-id`, sample_type, total_asvs, pct_flagged) %>%
  rename(
    "Sample ID"   = `sample-id`,
    "Sample Type" = sample_type,
    "Total ASVs"  = total_asvs,
    "% Flagged"   = pct_flagged
  )

# ==========================
# RENDER PNG
# ==========================

header_theme <- ttheme_default(
  colhead = list(
    bg_params = list(fill = "#2E86AB", col = "white"),
    fg_params = list(col = "white", fontface = "bold")
  ),
  core = list(
    bg_params = list(fill = c("white", "#F0F0F0")),
    fg_params = list(col = "black")
  )
)

warning_theme <- ttheme_default(
  colhead = list(
    bg_params = list(fill = "#E63946", col = "white"),
    fg_params = list(col = "white", fontface = "bold")
  ),
  core = list(
    bg_params = list(fill = c("white", "#FFF0F0")),
    fg_params = list(col = "black")
  )
)

grob_list <- list(
  textGrob("Contamination Summary Report",
           gp = gpar(fontsize = 16, fontface = "bold")),
  textGrob(paste0("Threshold for high contamination: ", CONTAMINATION_THRESHOLD, "%"),
           gp = gpar(fontsize = 11, col = "grey40")),
  textGrob("Group Averages (patient samples only)",
           gp = gpar(fontsize = 13, fontface = "bold")),
  tableGrob(group_averages, rows = NULL, theme = header_theme)
)

if (nrow(high_contam) > 0) {
  grob_list <- c(grob_list, list(
    textGrob(paste0("Samples above ", CONTAMINATION_THRESHOLD, "% contamination threshold"),
             gp = gpar(fontsize = 13, fontface = "bold", col = "black")),
    tableGrob(high_contam, rows = NULL, theme = warning_theme)
  ))
} else {
  grob_list <- c(grob_list, list(
    textGrob(paste0("No patient samples exceed the ", CONTAMINATION_THRESHOLD, "% contamination threshold."),
             gp = gpar(fontsize = 12, col = "#2E86AB", fontface = "italic"))
  ))
}

base_height <- 550

if (nrow(high_contam) > 0) {
  table2_height <- nrow(high_contam) * 60 + 150
} else {
  table2_height <- 80
}

png("contamination_summary.png",
    width  = 1000,
    height = base_height + table2_height,
    res    = 120)
grid.arrange(grobs = grob_list, ncol = 1)
dev.off()