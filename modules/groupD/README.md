# Group D — Visualisation Module 

## What this module does

This section of the pipeline takes the processed and validated outpus from previous steps, see group_B and group_C readme files for more information, and consolidates the data on a patient and population level by generating reports and graphs for clinicians. The purpose of this output is ultimately to help clinicians easily identify how their patients differ from an average healthy population, and what taxa the 16s sequencing reveals. The plots also summarize the population level symptoms of the patients. 

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

## What the report shows (for biologists) --> UPDATE WITH SPECIFIC CHARTS AND HOW TO INTERPRET. 

### Taxonomy
**What it answers:** Which bacteria are most abundant, and does their composition differ between treatment groups?

- **Stacked bar chart** — shows relative abundance of each bacterial family for every sample. Use the filter buttons to compare T0 vs T84, or one treatment group vs another.
- **Top taxa ranking** — horizontal bar showing which families dominate on average.
- **Sunburst** — hierarchical view of group → family → abundance.

### Alpha Diversity
**What it answers:** How rich and even is the microbial community within each person?
- **Rarefaction curves** — shows whether sequencing depth was sufficient to capture community richness.
- **Multi-metric chart** — observed richness (bars) and Pielou J′ evenness (diamonds) side by side.

### Beta Diversity
**What it answers:** How different are the microbial communities between people or between timepoints?

- **PCoA (Bray-Curtis and Jaccard)** — ordination plots where closer dots mean more similar communities. If samples from different groups cluster separately, supplementation may have structured the microbiome.
- **NMDS** — alternative ordination emphasising rank-order distances.
- **Hierarchical dendrogram** — which samples are most similar to each other? Branches that cluster by group or timepoint indicate a treatment effect.

### Individual Analysis
**What it answers:** How did each patient's microbiome change individually?
- **Virulence Plot** - This graph looks at the relative abundance of each patient and highlights if a predetermined set of virulent bacteria have an abnormally high relative % in the patient compared to the control. 

### Comparative
**What it answers:** Which bacterial families significantly changed in abundance?
- **Volcano plot** — combines effect size (fold change) with statistical significance (FDR-corrected p-value). Red dots passed both thresholds.

---

## Running via Nextflow #### INPUT!!!!!



