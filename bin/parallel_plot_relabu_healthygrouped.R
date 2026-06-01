#!/usr/bin/env Rscript

# ============================================================
#  PARALLEL COORDINATES PLOT: Relative Abundance per Sample
#  - Healthy (yes): grouped in groups of ~GROUP_SIZE samples
#    (shows the group mean as a line)
#  - Healthy (no): all individual samples
#  Only samples with info in metadata
# ============================================================

library(tidyverse)
library(GGally)
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

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

PATH_ASV_TABLE <- opt$asv
PATH_TAXONOMY  <- opt$taxonomy
PATH_METADATA  <- opt$metadata

# --- Parameters ---------------------------------------------
TOP_N_GENERA   <- 15
GROUP_SIZE     <- 5        # Samples per group in the "yes" group
OUTPUT_FILE    <- "parallel_plot_relative_abundance_healthygrouped.pdf"
PLOT_WIDTH     <- 18
PLOT_HEIGHT    <- 8
COLOR_BY       <- "healthy"
ALPHA_LINES    <- 0.55
LINE_SIZE      <- 0.6

# --- 1. Load data ----------------------------------------
message("Loading data...")
asv_table <- read_tsv(PATH_ASV_TABLE, show_col_types = FALSE) %>%
  rename(ASV_ID = `#OTU ID`)
taxonomy <- read_tsv(PATH_TAXONOMY, show_col_types = FALSE) %>%
  select(ASV_ID, Genus, Family)
metadata <- read_tsv(PATH_METADATA, show_col_types = FALSE)

# --- 2. Valid samples (with healthy info) -----------------
id_col <- if ("sample-id" %in% colnames(metadata)) "sample-id" else "sample_name"

valid_meta <- metadata %>%
  filter(!is.na(healthy)) %>%
  select(all_of(c(id_col, COLOR_BY))) %>%
  rename(sample = all_of(id_col)) %>%
  filter(sample %in% colnames(asv_table)) %>%
  mutate(healthy = as.factor(healthy))

message(paste("Valid samples:", nrow(valid_meta),
              "| yes:", sum(valid_meta$healthy == "yes"),
              "| no:", sum(valid_meta$healthy == "no")))

# --- 3. Relative abundance at Genus level -------------------
message("Calculating relative abundances...")
asv_filt <- asv_table %>%
  select(ASV_ID, all_of(valid_meta$sample))

genus_rel <- asv_filt %>%
  pivot_longer(-ASV_ID, names_to = "sample", values_to = "counts") %>%
  left_join(taxonomy, by = "ASV_ID") %>%
  mutate(
    Genus = str_replace_all(coalesce(Genus, ""), "^g__", ""),
    Genus = if_else(Genus == "", paste0("Unclassified_", coalesce(Family, "Unknown")), Genus)
  ) %>%
  group_by(sample, Genus) %>%
  summarise(counts = sum(counts), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(rel_abund = counts / sum(counts)) %>%
  ungroup()

# --- 4. Top N genera ---------------------------------------
top_genera <- genus_rel %>%
  filter(Genus != "Unclassified_Unknown") %>%
  group_by(Genus) %>%
  summarise(mean_rel = mean(rel_abund), .groups = "drop") %>%
  arrange(desc(mean_rel)) %>%
  slice_head(n = TOP_N_GENERA) %>%
  pull(Genus)

# --- 5. Wide matrix (samples × genera) -------------------
plot_wide <- genus_rel %>%
  filter(Genus %in% top_genera) %>%
  select(sample, Genus, rel_abund) %>%
  pivot_wider(names_from = Genus, values_from = rel_abund, values_fill = 0) %>%
  select(sample, all_of(top_genera)) %>%
  inner_join(valid_meta, by = "sample")

# --- 6. Construction of the data.frame for plotting ------------
#   - "no" samples: all individual
#   - "yes" samples: grouped in blocks of GROUP_SIZE, using the mean

df_no <- plot_wide %>%
  filter(healthy == "no") %>%
  mutate(line_id = sample, healthy_label = "no (individual)")

# Assign groups to "yes" (random order to avoid biases)
set.seed(42)
yes_samples <- plot_wide %>%
  filter(healthy == "yes") %>%
  slice_sample(prop = 1)   # shuffle rows to avoid grouping by similar samples

n_yes    <- nrow(yes_samples)
n_groups <- ceiling(n_yes / GROUP_SIZE)

df_yes_grouped <- yes_samples %>%
  mutate(group_id = paste0("yes_g", ceiling(row_number() / GROUP_SIZE))) %>%
  group_by(group_id) %>%
  summarise(across(all_of(top_genera), mean), .groups = "drop") %>%
  mutate(
    line_id       = group_id,
    healthy_label = paste0("yes (mean ~", GROUP_SIZE, " samples)")
  )

message(paste("Lines 'no':", nrow(df_no),
              "| Groups 'yes':", nrow(df_yes_grouped)))

# Combine
df_plot <- bind_rows(
  df_no       %>% select(line_id, healthy_label, all_of(top_genera)),
  df_yes_grouped %>% select(line_id, healthy_label, all_of(top_genera))
) %>%
  mutate(healthy_label = factor(healthy_label,
           levels = c("no (individual)",
                      paste0("yes (mean ~", GROUP_SIZE, " samples)"))))

# --- 7. Column indices for numeric genera -----------
genus_col_idx <- which(colnames(df_plot) %in% top_genera)

# --- 8. Colores ---------------------------------------------
yes_label <- paste0("yes (mean ~", GROUP_SIZE, " samples)")
color_pal <- c(
  "no (individual)" = "#E41A1C"
)

color_pal[yes_label] <- "#377EB8"

# --- 9. Plot ------------------------------------------------
message("Generating parallel coordinates plot...")

p <- ggparcoord(
  data        = df_plot,
  columns     = genus_col_idx,
  groupColumn = "healthy_label",
  scale       = "uniminmax",
  alphaLines  = ALPHA_LINES,
  showPoints  = FALSE
) +
  geom_line(linewidth = LINE_SIZE) +
  scale_color_manual(values = color_pal, name = "healthy") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1, size = 10, face = "bold"),
    axis.text.y      = element_text(size = 8),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "grey40")
  ) +
  labs(
    title    = paste0("Parallel Coordinates: Relative Abundance (Top ", TOP_N_GENERA, " Genera)"),
    subtitle = paste0("No healthy: ", nrow(df_no), " individual lines | ",
                      "Healthy: ", nrow(df_yes_grouped), " groups (mean of ~",
                      GROUP_SIZE, " samples each)"),
    x = "Bacterial genus",
    y = "Normalized relative abundance"
  )

# PNG
ggsave(
  file.path(
    opt$outdir,
    "parallel_plot_relative_abundance_healthygrouped.png"
  ),
  plot = p,
  width = PLOT_WIDTH,
  height = PLOT_HEIGHT,
  dpi = 300,
  device = "png"
)

message("Parallel plot saved as PNG")
