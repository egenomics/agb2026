# Group D — Visualisation Module 

## What this module does

This section of the pipeline takes the processed and validated outpus from previous steps, see group_B and group_C readme files for more information, and consolidates the data on a patient and population level by generating reports and graphs for clinicians. The purpose of this output is ultimately to help clinicians easily identify how their patients differ from an average healthy population, and what taxa the 16s sequencing reveals. The plots also summarize the population level symptoms of the patients. 

Reports are split into two subdirectories, Exploratory charts, and Explanatory Charts. Please see those sections for more information. 

## I/O  -- The following tables show the data inputs and the following outputs of this section of the pipeline

### Input files
| Input | Description |
|:---:|:---:|
| `shannon/alpha-diversity.tsv` | Shannon entropy per sample — measures species diversity and evenness within each sample |
| `observed/alpha-diversity.tsv` | Observed features per sample — raw count of unique ASVs detected in each sample |
| `faith/alpha-diversity.tsv` | Faith's phylogenetic diversity per sample — accounts for evolutionary relatedness of detected ASVs |
| `simpson/alpha-diversity.tsv` | Simpson index per sample — measures dominance; values close to 1 indicate no single species dominates |
| `weighted_unifrac/distance-matrix.tsv` |  Pairwise dissimilarity between samples weighted by abundance — used to detect community-level differences |
| `bray_curtis/distance-matrix.tsv` |  Pairwise dissimilarity between samples based on abundance only — no phylogenetic tree required |
| `rarefaction_plots.png` | Shows whether sequencing depth was sufficient to capture the full diversity of each sample |
| `contamination_summary.png` |  Reports potential contamination detected in samples, typically derived from blank/negative controls |
| `patient-sample-information.tsv` | Anonymous patient information regarding health and symptoms |
| `asv_taxonomy.tsv` | ASV values categorized by taxonomy |
| `asv_table.tsv` | ASV counts |

### Output Charts
| Chart Name | Description |
|:----------:|:------------|
| `alpha_diversity_status` | Compares alpha diversity metrics (Shannon, observed features, Faith's PD, Simpson) across samples to assess within-sample diversity |
| `overview_table` | Summary table of key metrics and sample-level statistics across the cohort |
| `PCoA.nf` | Principal Coordinates Analysis plot showing between-sample (beta) diversity using UniFrac and Bray-Curtis distance matrices |
| `violin_plot` | Violin plot comparing distribution of virulent genus abundance between healthy and non-healthy patient groups |
| `virulence_patientX` | Per-patient stacked bar plots showing relative abundance of virulent vs non-virulent taxa against healthy control average |
| `clinical_association_map` | Co-occurrence heatmap showing correlation between microbial taxa across samples |
| `demographic_table` | Demographic summary table linking patient metadata to sample information |
| `heatmap_samples_bacteria` | Heatmap of bacterial relative abundance across all samples |
| `parallel_plot_relative_abundance` | Parallel coordinates plot comparing relative abundance profiles between healthy and non-healthy patients |
| `pca_individuals` | PCA plot of bacterial composition across samples |
| `beta_diversity_tree` | Phylogenetic tree annotated with beta diversity information |
| `volcano_plot` | Volcano plot highlighting differentially abundant taxa between healthy and non-healthy groups |


---

## Directory layout and file structure

```
results/
├── exploratory/
│
│   ├── clinical_association_map
│   │   └── clinical_association_map.png
│   │
│   ├── demographic_table
│   │   └── demographic_table.png
│   │
│   ├── heatmap_results
│   │   └── heatmap_samples_bacteria.png
│   │
│   ├── parallel_results
│   │   └── parallel_plot_relative_abundance_healthygrouped.png
│   │
│   ├── pca_results
│   │   ├── pca_samples_bacteria.rds
│   │   ├── genus_counts.rds
│   │   ├── pca_scree_plot.png
│   │   ├── pca_individuals.png
│   │   ├── pca_variables.png
│   │   ├── pca_biplot_healthy_vs_disease.png
│   │   ├── pca_pc1_density.png
│   │   └── pca_pc1_boxplot.png
│   │
│   ├── volcano_plot
│   │   └── volcano_plot.html
│   │
│   └── beta_diversity_tree
│       └── beta_diversity_tree.png
│
└── explanatory/
    │
    ├── alpha_patient
    │   ├── patient1_distribution.png
    │   ├── patient2_distribution.png
    │   ├── patient3_distribution.png
    │   └── ...
    │
    ├── alpha_zscore_status.rds
    │
    ├── PCoA_patient
    │   ├── PCA_highlight_patient1.png
    │   ├── PCA_highlight_patient2.png
    │   ├── PCA_highlight_patient3.png
    │   └── ...
    │
    ├── overview_sample
    │   ├── patient1_metadata_report.tsv
    │   ├── patient1_top_genus.tsv
    │   ├── patient2_metadata_report.tsv
    │   ├── patient2_top_genus.tsv
    │   └── ...
    │
    └── Virulence_analysis
        ├── virulence_patient1.png
        ├── virulence_patient2.png
        ├── virulence_patient3.png
        └── ...
```
---
## Exploratory Charts
---
- The exploratory charts are primarily meant to allow the physician to take a quick look at the overall patient population. These representations show demographic information, metadata reports, PCA resutls, relative abundance across all samples, volcano plots, and a tree plot of the differences between sample alpha diversity. These reports should be used to identify general trends in the sample, and can be useful. For more information on how alpha and beta diversity were calculated and why they were chosen as metrics for clinical diagnosis, please see the repository wiki [here](link to wiki).

## Explanatory Charts
- The explanatory charts are to help clinicians take a deep dive into single patient samples. Every chart in in the explanatory folder is produced per unhealthy patient and compared to the average of the healthy group. The explanatory 
---

## Running via Nextflow #### INPUT!!!!!
- Each of the .nf files calls the data output from the previous group's pipeline channel. For group C, the following channels are called:
-- ch_alpha
-- ch_beta
-- ch_rarefaction

From group B:
-- meta
-- asv_table
-- asv_taxonomy


### Nextflow Workflow
- All scripts used to calculate the metrics for the plots are found in the bin folder, while the .nf scripts used to manage the pipeline are found in the modules/exploratory and modules/explanatory folders. The NEXTFLOW MAIN FILE ultimately calls all of the scripts in the modules directories and runs SCRIPT to summarize them into a downloadable HTML/PDF report.



### Technical information
Languages: ![R](https://img.shields.io/github/r-package/v/:user/:repo) and ![Python](https://img.shields.io/badge/Python-3.10-3776AB?style=flat&logo=python&logoColor=white)
Pipeline Organization with: 







