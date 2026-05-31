#!/usr/bin/env Rscript

library(tidyverse)
library(decontam)

# =====================
# PARSE ARGUMENTS
# =====================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: decontam_validation.R <asv_table> <taxonomy> <metadata> <ground_truth>")
}

counts_file       <- args[1]
taxonomy_file     <- args[2]
metadata_file     <- args[3]
ground_truth_file <- args[4]

output_dir <- "results/validation_outputs"

# ===============
# LOAD INPUT DATA
# ===============

counts_df    <- read_tsv(counts_file, show_col_types = FALSE) %>%
  rename(ASV_ID = `#OTU ID`)
taxonomy_df  <- read_tsv(taxonomy_file, show_col_types = FALSE) %>%
  rename(ASV_ID = 1)
metadata_df  <- read_tsv(metadata_file, show_col_types = FALSE) %>%
  mutate(across(everything(), as.character))
ground_truth <- read_tsv(ground_truth_file, show_col_types = FALSE)

# Extract ground truth genera (first word of species name)
ground_truth <- ground_truth %>%
  mutate(genus = str_split(species, " ") %>% map_chr(1))

# Exclude Deinococcus radiodurans — cannot be amplified by primers
excluded_species <- "Deinococcus radiodurans"
excluded_genus   <- "Deinococcus"

gt_species <- ground_truth$species[ground_truth$species != excluded_species]
gt_genera  <- ground_truth$genus[ground_truth$genus != excluded_genus]

# ========================
# BUILD SAMPLE METADATA
# ========================

sample_metadata <- metadata_df %>%
  select(`sample-id`, healthy)

blank_ids <- sample_metadata$`sample-id`[sample_metadata$healthy == "blank"]
mock_ids  <- sample_metadata$`sample-id`[sample_metadata$healthy != "blank"]

cat(sprintf("Blank samples:       %d\n", length(blank_ids)))
cat(sprintf("Mock community samples: %d\n", length(mock_ids)))

# =======================
# RUN DECONTAM
# =======================

asv_matrix <- counts_df %>%
  column_to_rownames("ASV_ID") %>%
  as.matrix() %>%
  t()

sample_names <- rownames(asv_matrix)
is_blank     <- sample_names %in% blank_ids

contamination_results <- isContaminant(asv_matrix, neg = is_blank,
                                       method = "prevalence", threshold = 0.1)

flagged_asvs  <- rownames(contamination_results[contamination_results$contaminant == TRUE, ])
passed_asvs   <- rownames(contamination_results[contamination_results$contaminant == FALSE, ])

cat(sprintf("\nDecontam flagged %d ASVs as contaminants.\n", length(flagged_asvs)))
cat(sprintf("Decontam passed  %d ASVs.\n", length(passed_asvs)))

# =======================
# CROSS-REFERENCE WITH GROUND TRUTH
# =======================

# Clean genus names (handle hyphenated genera like Escherichia-Shigella)
taxonomy_df <- taxonomy_df %>%
  mutate(Genus_clean = str_split(Genus, "-") %>% map_chr(1))

# For each flagged ASV, check if its genus matches a ground truth genus
flagged_taxonomy <- taxonomy_df %>%
  filter(ASV_ID %in% flagged_asvs) %>%
  mutate(
    in_ground_truth = Genus_clean %in% gt_genera,
    verdict = if_else(in_ground_truth,
                      "FALSE POSITIVE — ground truth genus incorrectly flagged",
                      "Correctly flagged — not in ground truth")
  )

passed_taxonomy <- taxonomy_df %>%
  filter(ASV_ID %in% passed_asvs) %>%
  mutate(
    in_ground_truth = Genus_clean %in% gt_genera,
    verdict = if_else(in_ground_truth,
                      "Correctly passed — ground truth genus",
                      "Passed — not in ground truth")
  )

# =======================
# SUMMARY STATISTICS
# =======================

n_false_positives <- sum(flagged_taxonomy$in_ground_truth)
n_correct_flags   <- sum(!flagged_taxonomy$in_ground_truth)
n_correct_passes  <- sum(passed_taxonomy$in_ground_truth)
n_total_gt_genera <- length(gt_genera)

cat(sprintf("\n--- Decontam Validation Summary ---\n"))
cat(sprintf("Ground truth genera (excl. Deinococcus): %d\n", n_total_gt_genera))
cat(sprintf("Ground truth genera correctly passed:    %d\n", n_correct_passes))
cat(sprintf("Ground truth genera incorrectly flagged: %d (FALSE POSITIVES)\n", n_false_positives))
cat(sprintf("Non-ground-truth ASVs correctly flagged: %d\n", n_correct_flags))

# =======================
# BUILD OUTPUT TABLE
# =======================

all_results <- bind_rows(
  flagged_taxonomy %>%
    select(ASV_ID, Genus, Genus_clean, in_ground_truth, verdict) %>%
    mutate(decontam_result = "Flagged"),
  passed_taxonomy %>%
    select(ASV_ID, Genus, Genus_clean, in_ground_truth, verdict) %>%
    mutate(decontam_result = "Passed")
) %>%
  arrange(decontam_result, desc(in_ground_truth), Genus_clean)

# Append summary rows
summary_rows <- tibble(
  ASV_ID          = c("SUMMARY", "SUMMARY", "SUMMARY", "SUMMARY"),
  Genus           = NA,
  Genus_clean     = NA,
  in_ground_truth = NA,
  decontam_result = NA,
  verdict         = c(
    paste0("Ground truth genera (excl. Deinococcus): ", n_total_gt_genera),
    paste0("Ground truth genera correctly passed: ",    n_correct_passes),
    paste0("Ground truth genera incorrectly flagged (FALSE POSITIVES): ", n_false_positives),
    paste0("Non-ground-truth ASVs correctly flagged: ", n_correct_flags)
  )
)

final_output <- bind_rows(all_results, summary_rows)

write_tsv(final_output, file.path(output_dir, "decontam_validation.tsv"))