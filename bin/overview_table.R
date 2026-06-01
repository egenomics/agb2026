library(tidyverse)
library(optparse)

arg_list <- list(
  # Metadata
  make_option(c("-m", "--metadata"), type = "character", help = "Path to sample-metadata.tsv", metavar = "path"),
  
  # RDS objects
  make_option(c("--pca"), type = "character", help = "Path to faith.tsv", metavar = "path"),
  make_option(c("--zscored"), type = "character", help = "Path to observed.tsv", metavar = "path"),
  make_option(c("--counts"), type = "character", help = "Path to shannon.tsv", metavar = "path"),
  
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

# Checking the existance of the files
files_to_check <- c("pca", "zscored", "counts", "outdir")

for (f in files_to_check) {
  path <- opt[[f]]
  if (!is.null(path)) {
    cat(sprintf("Loaded %s: %s\n", f, path))
  } else {
    cat(sprintf("Warning: %s was not provided.\n", f))
  }
}

# Create exist.ok
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

meta <- read_tsv(opt$metadata)
PrCA <- readRDS(opt$pca)
Zscores <- readRDS(opt$zscored)
genus_counts <- readRDS(opt$counts)


meta_filtered <- meta %>% filter(healthy == "no")
meta_filtered$healthy <- factor(meta_filtered$healthy, 
                                levels = c("no", "yes"), 
                                labels = c("Non-Healthy", "Healthy"))

# Genus counts only Non-healthy
genus_counts <- genus_counts %>%
  filter(sample %in% (meta_filtered %>% pull(`sample-id`)))

# 1. Extract the overall loadings for PC1 from the rotation matrix
pc1_overall_loadings <- as.data.frame(PrCA$rotation) %>%
  rownames_to_column("Feature") %>%
  select(Feature, PC1) %>%
  # 2. Calculate the absolute weight to find the strongest drivers
  mutate(Abs_Weight = abs(PC1)) %>%
  # 3. Grab the top 5 highest absolute weights
  slice_max(order_by = Abs_Weight, n = 10) %>%
  # 4. Clean up and sort
  arrange(desc(Abs_Weight))

# Adding the probability of this alpha Diversity
meta_filtered <- meta_filtered %>% mutate(
  alpha_diversity_probability = Zscores %>% pull(tail_prob)
)

# Covert every column into row
for (pid in meta_filtered %>% pull(`sample-id`)){
  
  sample_genus_counts <- genus_counts %>%
    filter(sample == pid) %>%
    arrange(counts) %>%
    arrange(-row_number())
  
  n_reads <- sample_genus_counts %>% pull(counts) %>% sum()
    
  sample_genus_counts <- sample_genus_counts %>% head(n = 5) %>% 
    mutate(
      relative_abundance = paste(round((counts/n_reads)* 100, digits = 3), "%")
    ) %>% mutate(
      in_Drivers = Genus %in% (pc1_overall_loadings %>% pull(Feature))
    )
  
  patient_metadata_long <- meta_filtered %>%
    filter(`sample-id` == pid) %>%
    # Convert all columns to character so they can coexist in one 'Value' column
    mutate(across(everything(), as.character)) %>%
    pivot_longer(
      cols = everything(), 
      names_to = "Attribute", 
      values_to = "Value"
    )

  metadata_path <- file.path(opt$outdir, paste0(pid, "_metadata_report.tsv"))
  genus_path    <- file.path(opt$outdir, paste0(pid, "_top_genus.tsv"))
  
  write_tsv(patient_metadata_long, metadata_path)
  write_tsv(sample_genus_counts, genus_path)
}
