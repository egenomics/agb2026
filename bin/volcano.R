#!/usr/bin/env Rscript

# 1. Fetch Nextflow Arguments
args <- commandArgs(trailingOnly = TRUE)
asv_file  <- args[1]
tax_file  <- args[2]
meta_file <- args[3]
out_dir   <- args[4] 

library(ANCOMBC)
library(phyloseq)
library(tidyverse)

# Read inputs directly from arguments
asv_df <- read.delim(asv_file, row.names = 1, sep = "\t", check.names = FALSE)
tax_df <- read.delim(tax_file, sep = "\t", check.names = FALSE)
meta_df <- read.delim(meta_file, sep = "\t", check.names = FALSE)

# 2. DATA PROCESSING & PHYLOSEQ OBJECT CREATION
meta_df <- meta_df %>%
  mutate(healthy_status = ifelse(healthy == "yes", "Healthy", "Has Condition")) %>%
  mutate(healthy_status = factor(healthy_status, levels = c("Healthy", "Has Condition")))

asv_samples <- colnames(asv_df)


meta_df_valid <- meta_df %>% filter(`sample-id` %in% asv_samples)
overlap_samples <- meta_df_valid$`sample-id`

asv_subset <- asv_df[, overlap_samples, drop = FALSE]
asv_subset <- asv_subset[rowSums(asv_subset) > 0, , drop = FALSE]

rownames(meta_df_valid) <- meta_df_valid$`sample-id`
rownames(tax_df)        <- tax_df$ASV_ID

tax_mat <- as.matrix(tax_df[rownames(asv_subset), c("Family", "Genus")])

otu_obj  <- otu_table(as.matrix(asv_subset), taxa_are_rows = TRUE)
tax_obj  <- tax_table(tax_mat)
meta_obj <- sample_data(meta_df_valid)

ps <- phyloseq(otu_obj, tax_obj, meta_obj)

# 3. RUN ANCOM-BC2 DIFFERENTIAL ABUNDANCE
ancom_output <- ancombc2(
  data = ps, tax_level = NULL, fix_formula = "healthy_status",
  group = "healthy_status", p_adj_method = "BH", alpha = 0.05,
  struc_zero = TRUE, neg_lb = TRUE, verbose = FALSE
)

res_df <- ancom_output$res

# 4. PREPARE VOLCANO PLOT DATA FRAME
lfc_col <- "lfc_healthy_statusHas Condition"
p_col   <- "p_healthy_statusHas Condition"
q_col   <- "q_healthy_statusHas Condition"

plot_data <- res_df %>%
  transmute(ASV_ID = taxon, log2FoldChange = .data[[lfc_col]], pvalue = .data[[p_col]],
            qvalue = .data[[q_col]], negLog10P = -log10(pvalue)) %>%
  left_join(tax_df, by = c("ASV_ID" = "ASV_ID")) %>%
  mutate(
    Significance = case_when(
      qvalue < 0.05 & log2FoldChange > 1  ~ "Enriched in Condition",
      qvalue < 0.05 & log2FoldChange < -1 ~ "Depleted in Condition",
      TRUE                                ~ "Not Significant"
    )
  )

# 5. GENERATE AND SAVE STATIC PLOT
color_map <- c("Not Significant" = "#e0e0e0", "Enriched in Condition" = "#d9534f", "Depleted in Condition" = "#2b579a")

# Filter for ALL significant points (both Enriched and Depleted) for labeling
significant_points <- plot_data %>% filter(Significance != "Not Significant")

gg_volcano <- ggplot(plot_data, aes(x = log2FoldChange, y = negLog10P, color = Significance)) +
  geom_point(alpha = 0.7, size = 1.5) + 
  scale_color_manual(values = color_map) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.4) +
  geom_text(data = significant_points,
            aes(label = paste0(ASV_ID, " (", round(log2FoldChange, 1), ")")),
            size = 3, vjust = -1, check_overlap = TRUE, color = "black") +
  labs(title = "Differential Abundance (Volcano Plot)", x = "Log2 Fold Change", y = "-Log10 (P-Value)") +
  theme_minimal() + 
  theme(panel.background = element_blank(), panel.grid.major = element_line(color = "#f5f5f5"), panel.grid.minor = element_blank())


dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "volcano_plot.png"), plot = gg_volcano, width = 10, height = 8, dpi = 300)