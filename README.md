# Group C — System Validation & Diversity Analysis

## Overview

Group C ensures that the microbiome analysis system produces scientifically reliable,
reproducible, and robust results before downstream clinical interpretation and visualization.

Our work is divided into two major stages:

1. [System validation](#1-system-validation)
2. [Diversity Analysis & Depth Optimization](#2-diversity-analysis--depth-optimization)

---

## 1. System Validation

### Objective

Evaluate the performance and reliability of the microbiome analysis pipeline developed by
Group B using datasets with known microbial composition. This stage ensures that the pipeline
can accurately identify microbial taxa, estimate abundances, and remain stable under different
analytical conditions.

### Validation Datasets

#### Primary Dataset — BEI Mock Communities (PRJEB10949)

**Source:** Lluch et al. 2015, *PLoS ONE*  
**DOI:** https://doi.org/10.1371/journal.pone.0142334  
**ENA accession:** [PRJEB10949](https://www.ebi.ac.uk/ena/browser/view/PRJEB10949)

Illumina MiSeq 16S V3-V4 paired-end sequencing of two BEI Resources mock communities
(even and staggered) and water-only negative controls. Selected on instructor recommendation.

**Why this dataset:**
- Known ground truth composition (20 bacterial species, concentrations documented in
  Supplementary Table S1 of the paper)
- Even community tests baseline detection accuracy (all species at 5%)
- Staggered community tests performance under realistic abundance imbalance (0.03%–27.3%)
- Negative controls (H2O blanks) enable contamination filtering validation

**Runs selected:**

| Run accession | Sample type | Read count | Purpose |
|---|---|---|---|
| ERR1049996 | BEI even mock — rep 1 | ~145,000 | Benchmarking |
| ERR1049997 | BEI even mock — rep 2 | ~151,000 | Benchmarking |
| ERR1049998 | BEI even mock — rep 3 | ~167,000 | Benchmarking |
| ERR1049999 | BEI staggered mock — rep 1 | ~163,000 | Adversarial test |
| ERR1050000 | BEI staggered mock — rep 2 | ~157,000 | Adversarial test |
| ERR1050001 | BEI staggered mock — rep 3 | ~146,000 | Adversarial test |
| ERR1049992 | H2O negative control 1 | ~156,000 | Contamination filtering |
| ERR1049993 | H2O negative control 2 | ~169,000 | Contamination filtering |
| ERR1049994 | H2O negative control 3 | ~153,000 | Contamination filtering |
| ERR1049995 | H2O negative control 4 | ~127,000 | Contamination filtering |
| ERR1049938 | H2O negative control 5 | ~8,500 | Contamination + low-depth stress test |
| ERR1049939 | H2O negative control 6 | ~17,000 | Contamination + low-depth stress test |
| ERR1049940 | H2O negative control 7 | ~17,000 | Contamination + low-depth stress test |

**Ground truth file:** `data/ground_truth/ground_truth_PRJEB10949.tsv`

> **Known limitation:** *Deinococcus radiodurans* cannot be amplified by the primers
> used in this study. A result of 0% for this species is expected and does not indicate
> pipeline failure (confirmed in the original paper).

---
### Validation Workflow

#### Step 1 — Pipeline Execution

Validation datasets are processed using the complete microbiome analysis pipeline developed
by Group B.

#### Step 2 — Taxonomic Validation

Pipeline outputs are compared against the expected microbial composition. Validation includes:

- Presence/absence of taxa
- Relative abundance estimation (optionally)
- Taxonomic classification accuracy

#### Step 3 — Performance Metrics

Quantitative metrics are computed to evaluate pipeline performance across two categories:

**Taxonomic Detection**

| Metric | Description |
|---|---|
| Accuracy | Overall classification correctness |
| Precision | True positive rate among predicted positives |
| Recall | True positive rate among actual positives |
| F1-score | Harmonic mean of precision and recall |

**Abundance Estimation**

| Metric | Description |
|---|---|
| Bray-Curtis dissimilarity | Compositional distance between profiles |
| RMSE | Root Mean Square Error of abundance estimates |

### Contamination Filtering

To reduce false-positive detections caused by laboratory or reagent contamination
("kit-ome"), filtering procedures are applied using:

- Negative controls
- Blank samples
- Known contaminant databases

This step improves the biological reliability of all downstream analyses.

### Stress Testing

Pipeline robustness is evaluated under varying conditions, including:

- Reduced sequencing depth
- Increased noise levels
- Read subsampling
- Variable read quality
- Artificial contamination scenarios

The goal is to assess the stability and reproducibility of taxonomic and abundance outputs
across adverse conditions.

### Validation Thresholds

The pipeline is considered validated only when predefined quality thresholds are met, including:

- High precision and recall for taxonomic detection
- Stable abundance estimation across replicates
- Low dissimilarity between expected and predicted profiles
- Robustness across subsampling conditions

---

## 2. Diversity Analysis & Depth Optimization

### Objective

Once the pipeline has been validated, Group C performs downstream  diversity analyses
on the processed outputs from Group B. This stage focuses on ensuring statistically fair and
biologically meaningful comparisons between samples.

### Input Data

| Description | Input |
|---|---|
| ASV/OTU abundance tables | table_counts.tsv |
| Taxonomic profiles | rep_seqs.fasta |
| Sample metadata | sample-metadata.tsv |

### Alpha Diversity Analysis

Within-sample diversity is measured using the following metrics:

- Shannon Diversity Index
- Simpson Index
- Observed Features
- Faith

### Beta Diversity Analysis

Between-sample community differences are evaluated using:

- Bray-Curtis dissimilarity
- UniFrac distances *(if phylogenetic information is available)*

### Rarefaction & Subsampling

Rarefaction analyses are conducted to:

- Evaluate sequencing depth sufficiency
- Identify optimal subsampling thresholds
- Ensure fair comparisons across samples with uneven sequencing depth

---

## Outputs

All deliverables generated by Group C are passed to **Group D** for visualization and
clinical interpretation.

| Deliverable | Format |
|---|---|
| Diversity metrics | Tabular (`.tsv`) |
| Rarefaction curves | Figures |
| Distance matrices | `.tsv` |

<img width="1920" height="1080" alt="Flowchart_groupC" src="https://github.com/user-attachments/assets/f0ae4181-c42b-4acc-b561-e0ac62939d29" />
