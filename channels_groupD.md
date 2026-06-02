# Channels for the Report Module

For the report module, it is interesting to know there are two workflows working in different ways:
- Exploratory Module: One instance for all the samples
- Explanatory Module: One instance for each Non-Healthy patient

Each of these workflows has different outputs, those are the following

> How are outputs specified in Nexflow?
- The emitted output of a process or workflow can be accesed using the .out. method
- The name specified below are the same for the outputs emitted

## Exploratory Report

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

## Explanatory Report

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
