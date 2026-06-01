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

## Directory layout and file structure of outputs 

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
# Channels for the Report Module

For the report module, it is interesting to know there are two workflows working in different ways:
- Exploratory Module: One instance for all the samples
- Explanatory Module: One instance for each Non-Healthy patient

Each of these workflows has different outputs, those are the following

> How are outputs specified in Nexflow?
- The emitted output of a process or workflow can be accesed using the .out. method
- The name specified below are the same for the outputs emitted

## Exploratory Report
- The exploratory charts are primarily meant to allow the physician to take a quick look at the overall patient population. These representations show demographic information, metadata reports, PCA resutls, relative abundance across all samples, volcano plots, and a tree plot of the differences between sample alpha diversity. These reports should be used to identify general trends in the sample, and can be useful. For more information on how alpha and beta diversity were calculated and why they were chosen as metrics for clinical diagnosis, please see the repository wiki [here](https://github.com/egenomics/agb2026/wiki/Output).

Inside `/exploratory` folder -> **channel_name** (`generated_directory/`):

- **demo_table_png** (`demographic_table/`): Demographics summary table of the study cohort.
  - *demographic_table.png*: Distribution of demographic and lifestyle variables stratified by health status.
- **co_heat** (`clinical_association_map/`): Clinical metadata overview.
  - *clinical_association_map.png*: Heatmap showing the co-occurrence of clinical conditions across samples.
- PCA results:
    - **pca_results** (`pca_results/`): Principal Component Analysis (PCA) outputs based on genus-level microbiome composition.
      - *pca_biplot_healthy_vs_disease.png*: PCA biplot showing sample clustering according to health status.
      - *pca_variables.png*: Most influential bacterial genera driving variation in the first principal components.
      - *pca_pc1_boxplot.png*: Comparison of PC1 score distributions between healthy and non-healthy samples.
    These 3 PNGs are generated but are not included in the final report:
      - *pca_individuals.png*: PCA projection of samples coloured by their contribution quality (cos²).
      - *pca_scree_plot.png*: Percentage of variance explained by each principal component.
      - *pca_pc1_density.png*: Distribution of sample scores along PC1 stratified by health status.
    - **pca_object** (`pca_results/`): PCA object exported for downstream explanatory analyses.
      - *pca_samples_bacteria.rds*: The RDS object of the PCA for its use at the explanatory report.
    - **pca_genus_counts** (`pca_results/`): Genus-level abundance table exported for downstream explanatory analyses.
      - *genus_counts.rds*: Aggregated genus count table used to generate the PCA.
- **parallel_results** (`parallel_results/`): Relative abundance overview of dominant bacterial genera.
  - *parallel_plot_relative_abundance_healthygrouped.png*: Parallel plot showing normalized relative abundances of the most abundant genera across samples, comparing individual non-healthy samples against grouped healthy controls.
- **heatmap_results** (`heatmap_results/`): Global microbiome composition heatmap.
  - *heatmap_samples_bacteria.png*: Hierarchical clustering heatmap of the most abundant bacterial genera across samples using relative abundances.
- **volcano_png** (`volcano_plot/`): Differential abundance analysis.
  - *volcano_plot*: Interactive volcano plot showing differential abundance results, highlighting taxa enriched or depleted in non-healthy samples.
- **tree_png** (`beta_diversity_tree/`): Beta-diversity hierarchical clustering analysis.
  - *beta_diversity_tree.png*: Hierarchical clustering tree (UPGMA) generated from Bray-Curtis beta-diversity distances, with samples coloured according to health status.

--- 

## Explanatory Report
- The explanatory charts are to help clinicians take a deep dive into single patient samples. Every chart in in the explanatory folder is produced per unhealthy patient and compared to the average of the healthy group. To find more information on how the eplanatory plots were created, please see the wiki [here]((https://github.com/egenomics/agb2026/wiki/Output). 
  
Inside `/explanatory` folder -> **channel_name** (`generated_directory/`):

- **overview_sample_table_dir** (`overview_sample/`): Folder containing all the overview tables for each Non-Healthy patient.
  - *{patient_id}_metadata_report.tsv*: Metadata information for the selected patient in long format.
  - *{patient_id}_top_genus.tsv*: Top 5 most abundant bacterial genera detected in the selected patient.
- **alpha_div_dist_dir** (`alpha_patient/`): Folder containing all the alpha distributions for each Non-Healthy patient.
  - *{patient_id}_distribution.png*: Position of the selected non-healthy patient within the distribution of healthy samples based on standardized alpha-diversity metrics.
- **zscore_status**: Alpha-diversity statistics exported for downstream explanatory analyses.
  - *alpha_zscore_status.rds*: Table containing mean alpha-diversity Z-scores and associated tail probabilities for non-healthy patients.
- **PCoA_patient_plots_dir** (`PCoA_patient/`): Folder location of the highlight for each non-Healthy patient.
  - *PCA_highlight_{patient_id}.png*: PCA plot highlighting the position of a specific non-healthy patient relative to the healthy cohort.
- **virulence_report** (`Virulence_analysis/`): Virulence-associated taxonomic abundance analysis.
  - *virulence_{patient_id}.png*: Relative abundance comparison between a non-healthy patient and the average healthy control profile, highlighting predefined potentially virulent genera.
- **violin_report** (`violin_analysis/`): Patient-level abundance distribution visualizations.
  - *virulence_violin.png*: Violin plots showing the relative abundance distributions of detected virulent genera in healthy versus non-healthy samples, including Mann–Whitney U test significance statistics.

---
## Running via Nextflow
- All scripts used to calculate the metrics for the plots are found in the bin folder, while the .nf scripts used to manage the pipeline are found in the modules/exploratory and modules/explanatory folders. The `workflows/groupD.nf` ultimately calls all of the scripts in the modules directories and runs all scripts to summarize them into a downloadable HTML/PDF report. In the broader pipeline, this script is called by main.nf. Please see the root readme for information on how to run the full pipeline. 


Each of the .nf files calls the data output from the previous group's pipeline channel. From group_C, the following channels are called:
- ch_alpha
- ch_beta
- ch_rarefaction

From group_B:
- meta
- asv_table
- asv_taxonomy



### Technical information
Languages: 
- ![R](![GitHub R package version](https://img.shields.io/github/r-package/v/:user/:repo) 
- ![Python](https://img.shields.io/badge/Python-3.10-3776AB?style=flat&logo=python&logoColor=white)







