#!/usr/bin/Rscript

library(tidyverse)
library(optparse)

# Only takes the argument you pass

arg_list <- list(
  # Alpha Diversity (Diversity Table) Files
  make_option(c("--faith"), type = "character", help = "Path to faith.tsv", metavar = "path"),
  make_option(c("--observed"), type = "character", help = "Path to observed.tsv", metavar = "path"),
  make_option(c("--shannon"), type = "character", help = "Path to shannon.tsv", metavar = "path"),
  make_option(c("--simpson"), type = "character", help = "Path to simpson.tsv", metavar = "path"),
  
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

# Checking the existance of the files
files_to_check <- c("faith", "observed", "shannon", "simpson", "metadata", "outdir")

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

# Borrar after finishing pipeline
alpha_list =c(
  faith = opt$faith, 
  observed = opt$observed,
  shannon = opt$shannon,
  simpson = opt$simpson
  )

# Load metadata to separate by groups
meta <- invisible(read_tsv(opt$metadata))

rename_metric <- function(path, name){
  invisible(read_tsv(path)) %>%
    rename( !!"Sample" := 1, !!name := 2)
}

# Load and normalize Alphas
if(length(alpha_list) > 0){
  alpha_div <- imap(
    alpha_list, # List that is going to be applied a function
    ~rename_metric(.x, .y) # Executes the function: .x is the value of the list, .y is the name of the value
  ) %>% 
    reduce(full_join, by = "Sample")
  alpha_div_scaled <- alpha_div %>%
    mutate(across(where(is.numeric),
                  ~as.vector(scale(.)))) %>%
    mutate(
      mean = rowMeans(pick(-Sample), na.rm = TRUE)
    )
    
  cat("Alpha diversity tables merged and normalized.\n")
}

# Keep samples that are present only in the metadata
alpha_div_scaled <- alpha_div_scaled %>% filter(Sample %in% (meta %>% pull(`sample-id`)))
# Separate to keep the Healthy group
alpha_div_scaled_healthy <- alpha_div_scaled %>%
  filter(
    Sample %in%
      (meta %>% filter(healthy == "yes") %>% pull(`sample-id`)) # ADJUST TO REAL DATA
  )

# Creating the base plot with the distribution before iterating
# This plot indicates the alphas
healthy_dist <- ggplot(alpha_div_scaled_healthy, aes(x = mean))+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), 
                color = "grey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = c(-1.96, 1.96), color = "darkgrey", linetype = "dotted", linewidth = 1) +
  annotate("text", x = 1.96 - 0.35, y = 0.45, label = expression(alpha == 0.05), color = "black", hjust = 0) +
  annotate("text", x = -1.96 + 0.35, y = 0.45, label = expression(alpha == 0.05), color = "black", hjust = 1) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  theme_minimal()


#Iterating over each sample (non-healthy)
non_healthy_patients <- alpha_div_scaled %>% 
  filter(
    Sample %in%
      (meta %>% filter(healthy == "no") %>% pull(`sample-id`))
  )

# Extract the Z-score_table
alpha_p_value <- non_healthy_patients %>%
  mutate(
    # Get the area of the smaller tail (the probability/p-value)
    tail_prob = pmin(pnorm(mean), 1 - pnorm(mean))
  ) %>%
  select(Sample, mean, tail_prob)

rds_filename <- file.path("alpha_zscore_status.rds")
saveRDS(alpha_p_value, file = rds_filename)

cat(paste("\nRDS table saved to:", rds_filename))
cat("\nProcess complete.\n")

patient_ids <- non_healthy_patients %>% pull(Sample)
patient_plots <- patient_ids %>% walk(
  function(id){
    p_value <- alpha_div_scaled %>%
      filter(Sample == id) %>%
      pull(mean)
    
    p <- healthy_dist +
      geom_vline(xintercept = p_value, color = "red")+
      annotate("text", x = p_value + 0.3, y = 0.5, label = id, color = "red")+
      labs(
        title = paste("Alpha diversity status for", id),
        x = "Mean Z-Score (Alpha Diversity)",
        y = "Density"
        )
    filename <- file.path(opt$outdir, paste0(id, "_distribution.png"))
    ggsave(filename = filename, plot = p, width = 6, height = 4, bg = "white")
})

cat("\nOverview Table generation complete.")
