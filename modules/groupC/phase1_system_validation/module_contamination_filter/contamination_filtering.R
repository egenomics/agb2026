#!/usr/bin/env Rscript

library(tidyverse)
library(decontam)

# =====================
# PARSE ARGUMENTS
# =====================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("Usage: Rscript contamination_filtering.R <asv_table> <taxonomy> <kraken2_output> <blanks_metadata> <patients_metadata>")
}

counts_file    <- args[1]
taxonomy_file  <- args[2]
kraken_file    <- args[3]
blanks_file    <- args[4]
patients_file  <- args[5]

# =====================
# PREPARE DISCARD LISTS
# =====================

alien_taxa <- c(
  "Homo sapiens",
  "Thermus aquaticus"
)

# ===============
# LOAD INPUT DATA
# ===============

counts_df <- read_tsv(counts_file, show_col_types = FALSE) %>%
  rename(ASV_ID = `#OTU ID`)
taxonomy_df       <- read_tsv(taxonomy_file, show_col_types = FALSE)
kraken_df         <- read_tsv(kraken_file,
                               col_names = c("Status", "ASV_ID", "Taxonomy", "Length", "LCA"),
                               show_col_types = FALSE)
blanks_metadata   <- read_tsv(blanks_file,   show_col_types = FALSE)
patients_metadata <- read_tsv(patients_file, show_col_types = FALSE)

taxonomy_df <- taxonomy_df %>% rename(ASV_ID = 1)

# ========================
# BUILD SAMPLE METADATA
# ========================

sample_metadata <- bind_rows(
  blanks_metadata   %>% mutate(across(everything(), as.character)),
  patients_metadata %>% mutate(across(everything(), as.character))
) %>%
  select(sra_id, healthy) %>%
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
is_blank     <- sample_names %in% blanks_metadata$sra_id

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
    )) %>%
    relocate(Contamination_Flag, .after = ASV_ID)
}

annotated_counts_df   <- apply_flags(counts_df)
annotated_taxonomy_df <- apply_flags(taxonomy_df)

# ==========================
# SUMMARY TABLE PER SAMPLE
# ==========================

counts_long <- annotated_counts_df %>%
  pivot_longer(-c(ASV_ID, Contamination_Flag),
               names_to  = "sra_id",
               values_to = "count") %>%
  filter(count > 0) %>%
  left_join(sample_metadata, by = "sra_id")

contamination_summary <- counts_long %>%
  group_by(sra_id, sample_type) %>%
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
# EXPORT
# ==========================

write_tsv(annotated_counts_df,   "annotated_table_counts.tsv")
write_tsv(annotated_taxonomy_df, "annotated_taxonomy.tsv")
write_tsv(contamination_summary, "contamination_summary.tsv")
