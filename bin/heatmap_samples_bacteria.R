# ============================================================
#  HEATMAP: Samples × Bacteria (Relative Abundance)
#  Taxonomic level: Genus
#  Only samples present in the metadata table (with healthy info)
# ============================================================

# --- Packages -----------------------------------------------
library(tidyverse)
library(pheatmap)
library(viridis)
library(optparse)

arg_list <- list(
  make_option(c("--asv"), type = "character",
              help = "ASV abundance table"),
  make_option(c("--taxonomy"), type = "character",
              help = "Taxonomy table"),
  make_option(c("-m", "--metadata"), type = "character",
              help = "Metadata table"),
  make_option(c("-o", "--outdir"), type = "character",
              default = "./parallel_results",
              help = "Output directory")
)

opt_parser <- OptionParser(option_list = arg_list)
opt <- parse_args(opt_parser)

required <- c("asv","taxonomy","metadata")

for(r in required){
  if(is.null(opt[[r]])){
    print_help(opt_parser)
    stop(paste("Missing argument:", r))
  }
}

dir.create(opt$outdir,
           recursive = TRUE,
           showWarnings = FALSE)

PATH_ASV_TABLE <- opt$asv
PATH_TAXONOMY  <- opt$taxonomy
PATH_METADATA  <- opt$metadata

# --- Parameters ---------------------------------------------
TOP_N_GENERA   <- 40
OUTPUT_FILE    <- "heatmap_samples_bacteria.pdf"
PLOT_WIDTH     <- 20
PLOT_HEIGHT    <- 14
ANNOTATION_COL <- "healthy"   # Variable of metadata to annotate columns (NULL to disable)

# --- 1. Load data ----------------------------------------
message("Loading data...")
asv_table <- read_tsv(PATH_ASV_TABLE, show_col_types = FALSE) %>%
  rename(ASV_ID = `#OTU ID`)

taxonomy <- read_tsv(PATH_TAXONOMY, show_col_types = FALSE) %>%
  select(ASV_ID, Genus, Family, Phylum)

metadata <- read_tsv(PATH_METADATA, show_col_types = FALSE)

# --- 2. Filter samples: only those in metadata with healthy info and present in ASV table ---
id_col <- if ("sample-id" %in% colnames(metadata)) "sample-id" else "sample_name"

valid_samples <- metadata %>%
  filter(!is.na(healthy)) %>%          # Only with healthy info
  pull(all_of(id_col)) %>%
  intersect(colnames(asv_table))        # Only those present in the ASV table

message(paste("Valid samples (with metadata and in asv_table):", length(valid_samples)))

# Filter the ASV table to the valid samples
asv_filt <- asv_table %>%
  select(ASV_ID, all_of(valid_samples))

# --- 3. Clean and aggregate at the Genus level -----------------
message("Calculating relative abundances by genus...")

asv_long <- asv_filt %>%
  pivot_longer(-ASV_ID, names_to = "sample", values_to = "counts") %>%
  left_join(taxonomy, by = "ASV_ID") %>%
  mutate(
    Genus = str_replace_all(coalesce(Genus, ""), "^g__", ""),
    Genus = if_else(Genus == "", paste0("Unclassified_", coalesce(Family, "Unknown")), Genus)
  )

genus_rel <- asv_long %>%
  group_by(sample, Genus) %>%
  summarise(counts = sum(counts), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(rel_abund = counts / sum(counts)) %>%
  ungroup()

# --- 4. Select top N genera by mean relative abundance across samples ------
top_genera <- genus_rel %>%
  filter(Genus != "Unclassified_Unknown") %>% 
  group_by(Genus) %>%
  summarise(mean_rel = mean(rel_abund), .groups = "drop") %>%
  arrange(desc(mean_rel)) %>%
  slice_head(n = TOP_N_GENERA) %>%
  pull(Genus)

message("Top genera selected:")
message(paste(top_genera, collapse = ", "))

# --- 5. Build matrix (genera × samples) ---------------
heatmap_mat <- genus_rel %>%
  filter(Genus %in% top_genera) %>%
  select(sample, Genus, rel_abund) %>%
  pivot_wider(names_from = sample, values_from = rel_abund, values_fill = 0) %>%
  column_to_rownames("Genus") %>%
  as.matrix()

# Ensure row order by mean abundance descending
heatmap_mat <- heatmap_mat[top_genera, ]

# Check for any NA, NaN, or Inf values in the heatmap matrix
stopifnot(all(is.finite(heatmap_mat)))

# log10 transformation (with pseudocount) to handle zeros and reduce skewness
heatmap_mat_log <- log10(heatmap_mat + 1e-5)

# --- 6. Column annotation (samples) --------------------
annotation_col_df <- NULL
if (!is.null(ANNOTATION_COL) && ANNOTATION_COL %in% colnames(metadata)) {
  annotation_col_df <- metadata %>%
    filter(!is.na(healthy)) %>%
    select(all_of(c(id_col, ANNOTATION_COL))) %>%
    filter(.data[[id_col]] %in% colnames(heatmap_mat)) %>%
    distinct(.data[[id_col]], .keep_all = TRUE) %>%
    rename(group = all_of(ANNOTATION_COL)) %>%
    mutate(
      group = case_when(
        group == "yes" ~ "healthy",
        group == "no"  ~ "non-healthy",
        TRUE ~ as.character(group)
      ),
      group = as.factor(group)
    ) %>%
    column_to_rownames(id_col)
}

annotation_colors <- list(
  group = c(
    "healthy"     = "#377EB8",  # azul
    "non-healthy" = "#E41A1C"   # rojo
  )
)

# --- 7. Generate heatmap --------------------------------------
message("Generating heatmap...")

heatmap_obj <- pheatmap(
  heatmap_mat_log,
  color             = viridis(100),
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  clustering_method = "ward.D2",
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  annotation_col    = annotation_col_df,
  annotation_colors = annotation_colors,
  fontsize_row      = 9,
  fontsize_col      = 5,
  angle_col         = 90,
  fontsize          = 8,
  main              = paste0(
    "Heatmap: Top ", TOP_N_GENERA,
    " genera most abundant\n(log10 relative abundance, n=",
    ncol(heatmap_mat), " samples)"
  ),
  border_color      = NA,
  cellheight        = 14,
  silent            = TRUE
)

# PNG
png(
  file.path(
    opt$outdir,
    "heatmap_samples_bacteria.png"
  ),
  width = PLOT_WIDTH,
  height = PLOT_HEIGHT,
  units = "in",
  res = 300
)

grid::grid.newpage()
grid::grid.draw(heatmap_obj$gtable)
dev.off()

message("Heatmap saved as PNG")
