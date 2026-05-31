# Group D — Visualisation Module (MicroSee Report Generator)

## What this module does

Takes QIIME2 TSV exports from upstream groups and generates a **single self-contained HTML report** (`microsee_report.html`). The report opens in any browser with no server, no installs, and **no internet required** — Plotly.js (4.3 MB) is embedded directly in the file, making it fully offline-compatible on HPC nodes.

## INPUTS and OUTPUTS GO HERE. MAKE A TABLE.

---

## Module layout

```
describe tree here
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



