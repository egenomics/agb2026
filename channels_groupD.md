# Channels for the Report Module

For the report module, it is interesting to know there are two workflows working in different ways:
- Exploratory Module: One instance for all the samples
- Explanatory Module: One instance for each Non-Healthy patient

Each of these workflows has different outputs, those are the following

> How are outputs specified in Nexflow?
- The emitted output of a process or workflow can be accesed using the .out. method
- The name specified below are the same for the outputs emitted

## Exploratory Report

Inside /exploratory folder

- **demo_table_png**: Demographics summary of the cohort
- **metadata_report**: Metadata report heatmap for all samples
- PCA results:
    - **pca_results**: The results provided by a PCA analysis in png format.
    - *pca_object*: RDS object for explanatory report.
    - *pca_genus_counts*: RDS object for the explanatory report.
- **parallel_plot**: plot of the relative abundance of bacterias across samples.
- **heatmap_plot**: Microbiome heatmap of all samples.
- **volcano_plot**: Microbiome volcano plot for all samples.
- **beta_tree**: Tree plot of the differences between sample alpha diversity.


## Explanatory Report

Inside /explanatory folder

- **overview_table**: Folder containing all the overview tables for each Non-Healthy patient.
    - *{patient_id}_metadata_report.tsv*: Contains the metadata of that specific patient in long format
    - *{patient_id}_top_genus.tsv*: Contains the top 5 genus represented in a patient
- **alpha_plots**: Folder containing all the alpha distributions for each Non-Healthy patient.
    - *{patient_id}_distribution.png*: Contains the location of a Non-Healthy patient in the distribution of healthy ones.
- pcoa_plots: Folder location of the highlight for each non-Healthy patient.
    - *{patient_id}_PCA_highlight.png*: Contains the highlight of the position of a Non-Healthy patient in the PCoA.
- **patogeny_plots**: 
- **violin_plots**: 