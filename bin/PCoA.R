#!/usr/bin/env Rscript

library(tidyverse)
library(optparse)
library(factoextra)
library(FactoMineR)
library(ggrepel)

# Only takes the argument you pass
arg_list <- list(
  # Beta Diversity Files
  make_option(c("--pca"), type = "character", help = "Path to the PCA class instance", metavar = "path"),
  # Metadata
  make_option(c("-m", "--metadata"), type = "character", help = "Path to sample-metadata.tsv", metavar = "path"),
  
  # Output directory
  make_option(c("-o", "--outdir"), type = "character", default = "./results", 
              help = "Output directory [default= %default]", metavar = "path")
)

# Parsing it
opt_parser <- OptionParser(option_list = arg_list)
opt <- parse_args(opt_parser)

# Validating that metadata exists
if (is.null(opt$metadata)) {
  print_help(opt_parser)
  stop("Missing metadata file! Use -m or --metadata.", call. = FALSE)
}
cat("--- Processing Files ---\n")

files_to_check <- c("pca", "metadata", "outdir")

for (f in files_to_check) {
  path <- opt[[f]]
  if (!is.null(path)) {
    cat(sprintf("Loaded %s: %s\n", f, path))
  } else {
    cat(sprintf("Warning: %s was not provided.\n", f))
  }
}

PrCA <- readRDS(file = opt$pca)

# Create exist.ok
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# Load metadata to separate by groups
meta <- invisible(read_tsv(opt$metadata)) %>%
  select(healthy, `sample-id`)

meta_filtered <- meta %>% filter(`sample-id` %in% rownames(PrCA$x))

meta_filtered$healthy <- factor(meta_filtered$healthy, 
                                levels = c("no", "yes"), 
                                labels = c("Non-Healthy", "Healthy"))

# The ID of the patient you want to label
non_h_ids <- meta_filtered %>% 
  filter(healthy == "Non-Healthy") %>% 
  pull(`sample-id`) %>%
  as.character()
rownames(PrCA$x) <- as.character(meta_filtered$`sample-id`)

# Checking how the PCA can depict the variance between groups
p <- fviz_pca_biplot(
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

for (pid in non_h_ids){
  label_data <- as.data.frame(PrCA$x[pid, , drop = FALSE])
  label_data$sample_id <- pid

  if (nrow(label_data) == 0) next

  p_pacient <- p +
    geom_point(
      data = label_data,
      aes(x = PC1, y = PC2),
      color = "red",       
      shape = 21,            
      size = 5,              
      stroke = 1.5           
    ) +
    geom_text_repel(
    data = label_data,
    aes(x = PC1, y = PC2, label = sample_id),
    color = "red",
    fontface = "bold",
    size = 4,
    box.padding = 0.5,
    point.padding = 0.5,
  ) +
    labs(title = paste("PCA Clustering:", pid),
         subtitle = paste("Status: ", meta_filtered %>%
           filter(`sample-id` == pid) %>%
           pull(healthy) )
         )

  output_file <- file.path(opt$outdir, paste0(pid, "_PCA_highlight", ".png"))
  
  ggsave(filename = output_file, 
         plot = p_pacient, 
         width = 10, height = 7, dpi = 300, bg = "white")
}
