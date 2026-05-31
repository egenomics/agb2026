# ============================================================
#  PCA: Samples × Bacterial Genera
# ============================================================

library(optparse)
library(tidyverse)
library(ggrepel)
library(compositions)
library(FactoMineR)
library(factoextra)

# ============================================================
# Arguments
# ============================================================

arg_list <- list(
  make_option(c("--asv"), type = "character", help = "ASV abundance table"),
  make_option(c("--taxonomy"), type = "character", help = "Taxonomy table"),
  make_option(c("-m", "--metadata"), type = "character", help = "Metadata table"),
  make_option(c("-o", "--outdir"), type = "character", default = "./exploratory_results", help = "Output directory")
)

opt_parser <- OptionParser(option_list = arg_list)
opt <- parse_args(opt_parser)

# ============================================================
# Validation
# ============================================================

required <- c("asv", "taxonomy", "metadata")

for (r in required) {
  if (is.null(opt[[r]])) {
    print_help(opt_parser)
    stop(paste("Missing argument:", r))
  }
}

dir.create(opt$outdir,
           recursive = TRUE,
           showWarnings = FALSE)

cat("ASV:", opt$asv, "\n")
cat("Taxonomy:", opt$taxonomy, "\n")
cat("Metadata:", opt$metadata, "\n")
cat("Output:", opt$outdir, "\n")

PATH_ASV_TABLE <- opt$asv
PATH_TAXONOMY  <- opt$taxonomy
PATH_METADATA  <- opt$metadata

# --- Parameters ---------------------------------------------
OUTPUT_FILE    <- "pca_samples_bacteria.pdf"
PLOT_WIDTH     <- 10
PLOT_HEIGHT    <- 8
COLOR_BY       <- "healthy"       # Variable to color points
SHAPE_BY       <- NULL            # Variable to set point shape (NULL = no shape)
LABEL_POINTS   <- FALSE           # TRUE to display sample IDs
PSEUDOCOUNT    <- 0.5             # For CLR: replaces zeros before log transformation

# --- 1. Load data ----------------------------------------
message("Loading data...")
asv_table <- read_tsv(PATH_ASV_TABLE, show_col_types = FALSE) %>%
  rename(ASV_ID = `#OTU ID`)
taxonomy <- read_tsv(PATH_TAXONOMY, show_col_types = FALSE) %>%
  select(ASV_ID, Genus, Family)
metadata <- read_tsv(PATH_METADATA, show_col_types = FALSE)

# --- 2. Valid samples -------------------------------------
id_col <- if ("sample-id" %in% colnames(metadata)) "sample-id" else "sample_name"

valid_meta <- metadata %>%
  filter(!is.na(healthy)) %>%
  select(all_of(c(id_col, COLOR_BY))) %>%
  rename(sample = all_of(id_col)) %>%
  filter(sample %in% colnames(asv_table)) %>%
  mutate(across(all_of(COLOR_BY), as.factor))

message(paste("Valid samples:", nrow(valid_meta)))

# --- 3. Relative abundance at Genus level ----------------------
message("Calculating relative abundances at Genus level...")
asv_filt <- asv_table %>%
  select(ASV_ID, all_of(valid_meta$sample))

genus_counts <- asv_filt %>%
  pivot_longer(-ASV_ID, names_to = "sample", values_to = "counts") %>%
  left_join(taxonomy, by = "ASV_ID") %>%
  mutate(
    Genus = str_replace_all(coalesce(Genus, ""), "^g__", ""),
    Genus = if_else(Genus == "", paste0("Unclassified_", coalesce(Family, "Unknown")), Genus)
  ) %>%
  group_by(sample, Genus) %>%
  summarise(counts = sum(counts), .groups = "drop")

# Export the genus counts for potential use in other analyses
saveRDS(genus_counts, file = file.path(opt$outdir, "genus_counts.rds"))

# Matrix: samples (rows) × genera (columns)
count_mat <- genus_counts %>%
  pivot_wider(names_from = Genus, values_from = counts, values_fill = 0) %>%
  column_to_rownames("sample") %>%
  as.matrix()

message(paste("Matrix dimensions:", nrow(count_mat), "samples ×", ncol(count_mat), "genera"))

# Remove genera that are zero across all samples
count_mat <- count_mat[, colSums(count_mat) > 0]

# --- 4. CLR Transformation ----------------------------------
# CLR: log(x_i / geometric_mean(x)) on each sample
# Requires adding pseudocount to handle zeros
message("Applying CLR transformation...")
count_pseudo <- count_mat + PSEUDOCOUNT

# clr() of {compositions} operates on rows (each row = a composition)
clr_mat <- as.matrix(clr(count_pseudo))
rownames(clr_mat) <- rownames(count_mat)
colnames(clr_mat) <- colnames(count_mat)

# --- 5. PCA -------------------------------------------------
message("Calculating PCA...")
PrCA <- prcomp(clr_mat, scale. = FALSE, center = TRUE)
# scale.=FALSE because CLR already homogenizes variance across components

# Save the PCA object as RDS for explanatory plots
saveRDS(PrCA, file.path(opt$outdir, "pca_samples_bacteria.rds"))

# Results for variables
PCAvar <- get_pca_var(PrCA)
# Dataframe of contributions
contribPCA <- as.data.frame(PCAvar$contrib)

# Obtain the top 15 contributors to PC1
top_15_pc1 <- contribPCA %>%
  rownames_to_column(var = "Variable") %>%
  arrange(desc(Dim.1)) %>%
  slice_head(n = 15) %>%
  pull(Variable)

message("\n============================================================")
message(" Top 15 variables that contribute the most to PC1:")
message("============================================================")
message(paste(1:15, top_15_pc1, sep = ". ", collapse = "\n"))
message("============================================================\n")

# Vector with the drivers of the first 2 PCs named as the protein
vargroups <- ifelse(PCAvar$contrib[, "Dim.1"] > PCAvar$contrib[, "Dim.2"],
                    "PC1 Driver",
                    "PC2 Driver")
# Importances~/Documents/Workspace/DMI-WorkSpace/assignment-2-dtmk10
summary(PrCA)

# Visualize Eigenvalues (Scree Plot)
p_eig <- fviz_eig(PrCA)

ggsave(file.path(opt$outdir, "pca_scree_plot.png"), plot = p_eig, width = 8, height = 6, dpi = 300)

# Graph of Individuals (Samples)
p_ind <- fviz_pca_ind(
  PrCA,
  col.ind = "cos2",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  geom = "point",
  title = "Samples - PCA"
)

ggsave(file.path(opt$outdir, "pca_individuals.png"), plot = p_ind, width = 10, height = 8, dpi = 300)

# Graph of Variables (Arrow Plot) with the first 10 contributors
p_var <- fviz_pca_var(
  PrCA,
  repel = TRUE,
  select.var = list(contrib = 15),
  axes = c(1,2),
  col.var = vargroups,
  palette = c("green", "blue")
)

ggsave(file.path(opt$outdir, "pca_variables.png"), plot = p_var, width = 10, height = 8, dpi = 300)

meta_filtered <- metadata %>% filter(`sample-id` %in% rownames(PrCA$x))

# Checking how the PCA can depict the variance between groups
p_biplot <- fviz_pca_biplot(
  PrCA,
  col.var = "red",
  col.ind = as.factor(meta_filtered$healthy),
  addEllipses = TRUE,
  repel = TRUE,
  palette = c("#FC4E07", "#149202"),
  label = "none",
  geom.var = "none",
  invisible = "var",
  pointsize = 2
)

ggsave(file.path(opt$outdir, "pca_biplot_healthy_vs_disease.png"), plot = p_biplot, width = 10, height = 8, dpi = 300)

pc1_data <- data.frame(
  PC1 = PrCA$x[, 1],                           # Grabs the first column of the PCA coordinates
  Healthy = as.factor(meta_filtered$healthy)   # Your aligned grouping variable
)
pc1_variance <- round(summary(PrCA)$importance[2, 1] * 100, 1)

p_density <- ggplot(pc1_data, aes(x = PC1, fill = Healthy)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("#FC4E07", "#00AFBB")) +
  theme_minimal() +
  labs(
    title = "Distribution of Samples along PC1",
    x = paste0("PC1 (", pc1_variance, "% explained variance)"),
    y = "Density",
    fill = "Health Status"
  ) +
  theme(legend.position = "top")

ggsave(file.path(opt$outdir, "pca_pc1_density.png"), plot = p_density, width = 8, height = 6, dpi = 300)

p_boxplot <- ggplot(pc1_data, aes(x = Healthy, y = PC1, fill = Healthy)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 2) +
  scale_fill_manual(values = c("#FC4E07", "#00AFBB")) +
  theme_minimal() +
  labs(
    title = "PC1 Score by Health Status",
    y = paste0("PC1 (", pc1_variance, "% explained variance)"),
    x = "Health Status"
  ) +
  theme(legend.position = "none")

ggsave(file.path(opt$outdir, "pca_pc1_boxplot.png"), plot = p_boxplot, width = 6, height = 6, dpi = 300)