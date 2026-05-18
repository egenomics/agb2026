# Group D — Visualisation Module (MicroSee Report Generator)

## What this module does

Takes QIIME2 TSV exports from upstream groups and generates a **single self-contained HTML report** (`microsee_report.html`). The report opens in any browser with no server, no installs, and **no internet required** — Plotly.js (4.3 MB) is embedded directly in the file, making it fully offline-compatible on HPC nodes.

> **Context:** MicroSee started as a web app (`microsee/` at the repo root). This report generator is the evolved, portable form of that idea — the same interactive charts, but crystallised into one file you can email, archive, or open anywhere.

---

## Module layout

```
modules/groupD/
├── README.md
└── microsee_report/
    ├── pyproject.toml              ← Python package (pip install -e .)
    ├── environment.yml             ← Conda environment definition
    ├── main.nf                     ← Nextflow process (MICROSEE_REPORT)
    ├── tests/
    │   ├── data/                   ← Fixture TSVs (4 patients × 2 timepoints, all metrics)
    │   │   ├── feature-table.tsv
    │   │   ├── taxonomy.tsv
    │   │   ├── metadata.tsv        ← includes sixmwt + il18 clinical columns
    │   │   └── alpha-diversity.tsv ← all 5 metrics: shannon, simpson, observed, pielou, faith_pd
    │   ├── test_parsers.py         ← Unit tests for all parsers + integrate()
    │   └── test_charts.py          ← Smoke tests for chart builders + HTML rendering
    └── report_generator/           ← Python report engine
        ├── __init__.py
        ├── generate_report.py      ← CLI entry point  (microsee-report command)
        ├── parsers.py              ← QIIME2 TSV parsers
        ├── models.py               ← Pydantic v2 data models
        ├── requirements.txt        ← Pinned dependencies (for reference)
        ├── Dockerfile              ← Container image definition
        └── charts/
            ├── config.py            ← Colour palette and Plotly layout defaults
            ├── utils.py             ← Shared colour helpers (group + taxon palette)
            ├── preprocessing.py     ← Shared row helpers (get_patient_timepoints, sorted_timepoints, …)
            ├── metrics.py           ← METRIC_LABELS, metric_value(), pielou_evenness()
            ├── stats_helpers.py     ← Pure stat functions (Wilcoxon, MW, Welch t, Spearman, BH-FDR)
            ├── distances.py         ← Bray-Curtis, Jaccard, PCoA, clustering
            ├── taxonomy.py          ← Stacked bar (27 filter combos), donut, sunburst
            ├── alpha.py             ← Strip/box/violin, brackets, rarefaction, multi-metric
            ├── beta.py              ← PCoA, NMDS, dendrogram, Δ abundance heatmap
            ├── individual.py        ← Slopegraph, stability, rank plot, radar, small multiples
            ├── comparative.py       ← LFC bar, volcano, ANCOM-style CLR, heatmap, correlation
            ├── clinical.py          ← Clinical slopegraphs, Shannon scatter, taxa×clinical heatmap
            ├── stats.py             ← Wilcoxon / MW tests, LME trajectory, PERMANOVA, summary table
            ├── insights.py          ← Dynamic text insights generated from chart payloads
            ├── orchestrator.py      ← compute_chart_data() + ReportConfig (section selection)
            ├── renderer.py          ← Fills HTML templates (cohort + per-patient reports)
            ├── template.html        ← HTML/CSS/JS cohort report shell (all interactive controls)
            ├── patient_template.html← Per-patient HTML report shell
            └── plotly.min.js        ← Bundled Plotly.js v2.35.2 (MUST be committed to git)
```

---

## Quickstart (works out of the box)

```bash
# 1. Install once (from repo root)
pip install -e "modules/groupD/microsee_report"

# 2. Generate a report using the bundled fixture data
microsee-report \
    --feature-table modules/groupD/microsee_report/tests/data/feature-table.tsv \
    --taxonomy      modules/groupD/microsee_report/tests/data/taxonomy.tsv \
    --metadata      modules/groupD/microsee_report/tests/data/metadata.tsv \
    --alpha         modules/groupD/microsee_report/tests/data/alpha-diversity.tsv \
    --output        /tmp/microsee_demo.html

# 3. Open in browser
open /tmp/microsee_demo.html          # macOS
xdg-open /tmp/microsee_demo.html      # Linux
```

---

## Usage with your own data

```bash
# --alpha is optional but strongly recommended
microsee-report \
    --feature-table path/to/feature-table.tsv \
    --taxonomy      path/to/taxonomy.tsv \
    --metadata      path/to/metadata.tsv \
    --alpha         path/to/alpha-diversity.tsv \
    --output        microsee_report.html

# Generate one HTML per patient (plus the combined report)
microsee-report ... --mode all --output microsee_report.html
# Produces: microsee_report.html + microsee_report_Pat1.html, etc.
```

**`--mode` options:**

| Mode | Output |
|---|---|
| `cohort` (default) | One combined report for all samples |
| `patient` | One HTML per patient |
| `all` | Combined report + one per patient |

---

## What the report shows (for biologists)

### Taxonomy
**What it answers:** Which bacteria are most abundant, and does their composition differ between treatment groups?

- **Stacked bar chart** — shows relative abundance of each bacterial family for every sample. Use the filter buttons to compare T0 vs T84, or one treatment group vs another.
- **Top taxa ranking** — horizontal bar showing which families dominate on average.
- **Donut chart** — average composition per group at a glance.
- **Sunburst** — hierarchical view of group → family → abundance.

### Alpha Diversity
**What it answers:** How rich and even is the microbial community within each person?

- **Strip / Box / Violin charts** — distributions of diversity metrics per group. Toggle between Shannon H′ (information content), Simpson 1−D (dominance), Pielou J′ (evenness), Observed taxa (richness), and Faith PD (phylogenetic diversity).
- **Significance brackets** — automatically computed Wilcoxon (paired, T0→T84) and Mann-Whitney (between groups at T84) p-values drawn directly on the box chart.
- **Rarefaction curves** — shows whether sequencing depth was sufficient to capture community richness.
- **Multi-metric chart** — observed richness (bars) and Pielou J′ evenness (diamonds) side by side.

### Beta Diversity
**What it answers:** How different are the microbial communities between people or between timepoints?

- **PCoA (Bray-Curtis and Jaccard)** — ordination plots where closer dots mean more similar communities. If samples from different groups cluster separately, supplementation may have structured the microbiome.
- **NMDS** — alternative ordination emphasising rank-order distances.
- **Hierarchical dendrogram** — which samples are most similar to each other? Branches that cluster by group or timepoint indicate a treatment effect.
- **Δ Abundance heatmap** — shows which families increased (red) or decreased (blue) between T0 and T84 for each patient. Sorted by the family that changed most.

### Individual Analysis
**What it answers:** How did each patient's microbiome change individually?

- **Paired slopegraph** — one line per patient showing Shannon H′ at T0 and T84. Lines going up = increased diversity.
- **Stability score** — Bray-Curtis dissimilarity for each patient (0 = identical T0 and T84 samples, 1 = completely different). Shorter bars = more stable microbiome.
- **Diversity rank** — all samples ranked from lowest to highest Shannon H′.
- **Patient radar** — spider/web chart showing group mean composition at T0 (filled) vs T84 (dashed). Each axis is one bacterial family.
- **NMDS trajectories** — arrows showing where each patient's community moved in ordination space (T0 → T84). Short arrows = stable, long arrows = large shift.
- **Small multiples** — individual stacked bar charts for every patient showing T0 and T84 side by side.

### Comparative
**What it answers:** Which bacterial families significantly changed in abundance?

- **Log fold change bar** — which families increased or decreased the most (T84 vs T0)?
- **Volcano plot** — combines effect size (fold change) with statistical significance (FDR-corrected p-value). Red dots passed both thresholds.
- **ANCOM-style CLR** — compositionally-unbiased differential abundance test (CLR-transformed paired Wilcoxon). More robust than simple fold-change for compositional data.
- **Abundance heatmap** — sample × family matrix, colour = relative abundance.
- **Taxon correlation matrix** — which families tend to co-occur (positive correlation) or compete (negative)?

### Clinical *(only shown if sixmwt / il18 columns are in metadata)*
**What it answers:** Does microbiome diversity correlate with physical function or inflammation?

- **6MWT and IL-18 slopegraphs** — individual patient trajectories for the 6-minute walk test and IL-18 cytokine, with group mean dashed line.
- **Shannon vs 6MWT / IL-18 scatter** — Pearson correlation between diversity and each clinical outcome, with regression line and r/p annotation.
- **Taxa × Clinical Spearman heatmap** — which families correlate with improvement in 6MWT or reduction in IL-18? Stars indicate statistical significance.

### Longitudinal
**What it answers:** How did diversity change over the study period on average?

- **Shannon over time** — group mean per timepoint connected by a line.
- **LME-style trajectory** — mean ± 95% CI band with individual patient lines underneath. Toggle the metric with the Alpha Diversity buttons. Wilcoxon p-value is annotated.

### Statistics
**What it answers:** Do the observed differences pass formal statistical tests?

- **PERMANOVA table** — tests which factor (supplementation group, timepoint, individual) explains the most variance in community composition. R² = fraction of variance explained.
- **Diversity summary table** — mean ± SD of all five alpha metrics per group.

---

## Running via Nextflow

```bash
# With conda (default)
nextflow run workflows/groupD.nf -profile conda \
    --feature_table path/to/feature-table.tsv \
    --taxonomy      path/to/taxonomy.tsv \
    --metadata      path/to/metadata.tsv \
    --alpha         path/to/alpha-diversity.tsv \
    --outdir        results/

# With bundled test fixtures
nextflow run workflows/groupD.nf -profile test,conda

# On a SLURM cluster
nextflow run workflows/groupD.nf -profile slurm,conda \
    --feature_table ...
```

> **Alpha diversity file** — QIIME2 exports one metric per file (e.g. `shannon_entropy.tsv`,
> `faith_pd.tsv`). Merge them into a single TSV before passing to `--alpha`:
> ```bash
> paste shannon_entropy.tsv <(cut -f2- faith_pd.tsv) \
>       <(cut -f2- observed_features.tsv) > alpha-diversity.tsv
> ```
> Without `--alpha`, Shannon and Simpson are estimated from family-level abundances,
> which underestimates diversity. All other metrics (Faith PD, Pielou, Observed) will be absent.

> **Docker / Singularity** — the container `ghcr.io/agb2026/microsee:latest` must be built
> and pushed before using `-profile docker` or `-profile singularity`. Build it with:
> ```bash
> docker build -t ghcr.io/agb2026/microsee:latest \
>     modules/groupD/microsee_report/report_generator/
> docker push ghcr.io/agb2026/microsee:latest
> ```
> Until then use `-profile conda` instead.

---

## Running tests

```bash
# Install with test extras
pip install -e "modules/groupD/microsee_report[dev]"

# Run the full test suite
pytest modules/groupD/microsee_report/tests/ -v
```

Tests cover inline string fixtures and the realistic 4-patient fixture TSV files.

---

## HPC compatibility checklist

- [x] **Offline** — Plotly.js bundled in the HTML; no CDN calls at render time
- [x] **No display required** — pure CLI, no GUI or X11 needed
- [x] **Shared filesystem safe** — reads from staged Nextflow work directory (NFS/Lustre compatible)
- [x] **Pure Python stats** — no R, no MATLAB, no scipy required
- [x] **Container support** — Dockerfile provided; use `-profile singularity` or `-profile docker`
- [x] **Conda support** — `environment.yml` provided with pinned deps

---

## Troubleshooting

### `plotly.min.js` is missing — charts are blank

The Plotly.js file must be committed to git so HPC nodes (no internet) can use it.

```bash
curl -fsSL https://cdn.plot.ly/plotly-2.35.2.min.js \
     -o modules/groupD/microsee_report/report_generator/charts/plotly.min.js

git add modules/groupD/microsee_report/report_generator/charts/plotly.min.js
git commit -m "Bundle Plotly.js v2.35.2 for offline HPC use"
```

The generator will attempt a one-time auto-download if the file is missing, but this fails on most HPC nodes. If you see `RuntimeWarning: charts/plotly.min.js not found`, run the commands above.

---

### `conda not found` or `conda: command not found` on HPC

Conda is often available via a module system on HPC clusters. Try:

```bash
module load anaconda3   # or: module load miniconda
conda activate microsee
```

Or use the Singularity profile instead:
```bash
nextflow run workflows/groupD.nf -profile singularity ...
```

---

### Sample IDs don't match between files

If you see a warning like `Sample X in feature-table but not in metadata`, the sample IDs are inconsistent. Common causes:

| Problem | Fix |
|---|---|
| Trailing whitespace in TSV | `sed -i 's/[[:space:]]*$//' metadata.tsv` |
| Mixed `_` vs `-` in IDs | Standardise to one separator across all files |
| `sample-id` vs `sample_id` header | The parser accepts both, but check for typos |
| QIIME2 added `#SampleID` prefix | The parser strips `#` from comment lines automatically |

The report will still generate for the intersecting samples and log a warning listing the mismatched IDs.

---

## Input file format reference

All inputs are QIIME2 TSV exports. Comment lines starting with `#` are skipped automatically.

| File | Required columns |
|---|---|
| `feature-table.tsv` | First col = feature/OTU ID; remaining cols = sample IDs with integer read counts |
| `taxonomy.tsv` | `Feature ID`, `Taxon` (semicolon-separated lineage, e.g. `d__Bacteria;p__Firmicutes;...`) |
| `metadata.tsv` | `sample-id`, plus any of: `subject`/`patient`, `group`/`treatment`, `timepoint`/`time` |
| `alpha-diversity.tsv` | `sample-id`, then any of: `shannon_entropy`, `simpson`, `observed_features`, `faith_pd`, `pielou_evenness` |

Column names are matched by regex so minor variations (`shannon` vs `shannon_entropy`) are handled automatically. Clinical columns `sixmwt` and `il18` are optional — if present, the Clinical section is included in the report.

---

## Workflow diagram

```mermaid
flowchart TD
    FT[feature-table.tsv]            --> R[microsee-report\nPython CLI]
    TX[taxonomy.tsv]                  --> R
    MD[metadata.tsv]                  --> R
    AL[alpha-diversity.tsv\noptional] --> R
    PJ[plotly.min.js\nbundled]        --> R
    R --> HTML[microsee_report.html\n~5 MB · self-contained]
    HTML --> USER[Open in any browser\nno server · no internet · no installs\nHPC compatible]
```
