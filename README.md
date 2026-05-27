Repository for the AGB 2026 common class project.  
**Paper:** [Short-Term Ingestion of Essential Amino Acid Based Nutritional Supplements or Whey Protein Improves the Physical Function of Older Adults Independently of Gut Microbiome](https://pubmed.ncbi.nlm.nih.gov/38426663/)
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

**Ground truth file:** `data/ground_truth_PRJEB10949.tsv`

> **Known limitation:** *Deinococcus radiodurans* cannot be amplified by the primers
> used in this study. A result of 0% for this species is expected and does not indicate
> pipeline failure (confirmed in the original paper).

> **Known limitation:** Group B's pipeline classifies ASVs to genus level only. Species-level
> metrics are therefore not computed. Genus names containing hyphens (e.g. `Escherichia-Shigella`)
> are truncated to the first component before matching against the ground truth.

For download instructions see the [project wiki](https://github.com/egenomics/agb2026/wiki/Pipeline-Validation-Datasets).

---

### Validation Workflow

#### Step 1 — Pipeline Execution

Validation datasets are processed using the complete microbiome analysis pipeline developed
by Group B.

#### Step 2 — Taxonomic Validation

Pipeline outputs are compared against the expected microbial composition. Validation includes:

- Presence/absence of taxa at genus level
- Relative abundance estimation at genus level
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

Metrics are computed per replicate and averaged per community type (even and staggered).

**Abundance Estimation**

| Metric | Description |
|---|---|
| Bray-Curtis dissimilarity | Compositional distance between observed and expected profiles |
| RMSE | Root Mean Square Error of abundance estimates |

#### Outputs

| File | Description |
|---|---|
| `detection_metrics.tsv` | Precision, recall, F1, accuracy per replicate and averaged per community type |
| `abundance_metrics.tsv` | RMSE and Bray-Curtis per replicate and averaged per community type |

#### Script

See `modules/groupC/phase1_system_validation/dataset_validation/validation_metrics.py`.

---

### Contamination Filtering

To reduce false-positive detections caused by laboratory or reagent contamination
("kit-ome"), filtering procedures are applied in two steps:

1. **Statistical filtering (decontam):** ASVs are flagged as contaminants if they are
   significantly more prevalent in negative controls than in real samples, using the
   prevalence method with a threshold of 0.1.

2. **Taxonomic filtering (Kraken2):** ASVs are flagged if they are classified as known
   biological contaminants, currently *Homo sapiens* and *Thermus aquaticus*, based on
   the Kraken2 database at `/data/upfagb/u269238/kraken2_db`.

> **Known limitation:** The Kraken2 database does not include Chloroplast, Mitochondria,
> or Halomonas sequences. Detection of these contaminants is therefore not possible with
> the current database.

ASVs are not removed but annotated with a `Contamination_Flag` column:

| Flag | Meaning |
|---|---|
| `Passed` | Not flagged by either method |
| `Flagged_Kitome` | Flagged by decontam only |
| `Flagged_Alien` | Flagged by Kraken2 only |
| `Flagged_Both` | Flagged by both methods |

#### Outputs

| File | Description |
|---|---|
| `annotated_table_counts.tsv` | ASV count table with contamination flags |
| `annotated_taxonomy.tsv` | Taxonomy table with contamination flags |
| `contamination_summary.tsv` | Per-sample summary of flagged ASV counts and percentages, classified by sample type (Blank / Healthy / Non-healthy) |

#### Modules

See `modules/groupC/phase1_system_validation/module_kraken2/` for the Kraken2 Nextflow module
and `modules/groupC/phase1_system_validation/module_contamination_filter/` for the
contamination filtering Nextflow module and R script.

---

### Stress Testing

Pipeline robustness is evaluated using targeted stress scenarios designed to test how the
system behaves under abnormal, adverse, or extreme input conditions. The goal is not only
to determine whether the pipeline completes successfully, but also whether it fails safely,
reports clear errors, and avoids generating misleading downstream results.

This first stress-testing stage uses the current Group B outputs as input:

- `asv_table.tsv`
- `ASV_taxonomy.tsv`
- `rep_seqs.fasta`

These files are used to prepare small derived test inputs at the ASV table, taxonomy, and
metadata level. This allows Group C modules to be tested before full end-to-end FASTQ-level
stress tests are run on the cluster.

The initial stress-test scenarios include:

| Test ID | Scenario | Purpose | Expected behaviour |
|---|---|---|---|
| `ST00` | Valid subset control | Confirm that a small valid input runs correctly | The module completes and generates expected outputs |
| `ST01` | Zero-count sample | Test behaviour when one sample has no reads | The sample is excluded or clearly flagged |
| `ST02` | Very low sequencing depth | Test robustness with samples downsampled to very few reads | Low-depth samples are flagged or removed during rarefaction/depth filtering |
| `ST03` | Invalid count table | Test behaviour when the ASV table contains a non-numeric count | The module fails early with a clear parsing error |
| `ST04` | Single-taxon dominance | Test biologically extreme but valid input | The module completes and reports very low diversity |
| `ST05` | Metadata mismatch | Test sample identifier inconsistencies between metadata and ASV table | The mismatch is detected before downstream analysis |
| `ST06` | Artificial contamination spike | Test whether control-enriched ASVs can be detected as potential contaminants | Spiked ASVs are flagged or reported as suspicious |

For each stress test, the following information will be recorded:

- input files used
- expected behaviour
- executed command or module
- exit status
- whether outputs were generated
- whether the error or warning message was clear
- whether any silent failure occurred
- final status: `PASS`, `FAIL`, or `WARNING`

A stress test can be considered successful even if the pipeline fails, as long as the failure
is expected, occurs early, and produces an interpretable error message. The main failure mode
to avoid is silent execution that produces apparently valid but biologically misleading outputs.

This section includes a documented stress-test folder with scenario definitions, small 
derived input files, scripts to regenerate the test inputs, and a results template.

### Validation Thresholds

The pipeline is considered validated only when predefined quality thresholds are met, including:

- High precision and recall for taxonomic detection
- Stable abundance estimation across replicates
- Low dissimilarity between expected and predicted profiles
- Robustness across subsampling conditions


## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/abgtemplate for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
