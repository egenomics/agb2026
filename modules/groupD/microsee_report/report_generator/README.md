# report_generator

Python package that powers the MicroSee self-contained HTML report. It reads QIIME2 TSV exports and writes a single offline HTML file (~5 MB) with 30+ interactive Plotly charts.

> For full documentation — Nextflow usage, HPC checklist, chart list, input formats — see [`modules/groupD/README.md`](../../README.md).

---

## Install

```bash
# From the repo root — installs the microsee-report CLI command
pip install -e "modules/groupD/microsee_report"

# With test dependencies
pip install -e "modules/groupD/microsee_report[dev]"
```

---

## Usage

```bash
microsee-report \
    --feature-table feature-table.tsv \
    --taxonomy      taxonomy.tsv      \
    --metadata      metadata.tsv      \
    [--alpha        alpha-diversity.tsv] \
    [--output       microsee_report.html]
```

`--alpha` is optional but unlocks Faith PD and more precise Pielou J′. Open the output file in any browser.

---

## Files

| File | Purpose |
|---|---|
| `generate_report.py` | CLI entry point (`microsee-report` command) |
| `parsers.py` | QIIME2 TSV parsing — feature-table, taxonomy, metadata, alpha-diversity |
| `integrator.py` | Joins all parsed outputs into chart-ready `SampleRow` objects |
| `models.py` | Pydantic v2 data models |
| `py.typed` | PEP 561 marker — declares this package as typed |
| `requirements.txt` | Pinned dependencies (reference; `pyproject.toml` is authoritative) |
| `charts/config.py` | Colour palette and Plotly layout/config defaults |
| `charts/utils.py` | Group + taxon colour helpers (cycling 20-colour palette for unknowns) |
| `charts/distances.py` | Bray-Curtis, Jaccard, PCoA, average-linkage clustering |
| `charts/taxonomy.py` | Stacked bar (27 filter variants), top-N ranking, donut, sunburst |
| `charts/alpha.py` | Strip/box/violin, significance brackets, rarefaction curves, multi-metric |
| `charts/beta.py` | PCoA (Bray-Curtis + Jaccard), NMDS, dendrogram, Δ abundance heatmap |
| `charts/individual.py` | Paired slopegraph, stability bar, rank plot, patient radar, small multiples |
| `charts/comparative.py` | LFC bar, volcano (BH-FDR), ANCOM-style CLR, abundance heatmap, correlation matrix |
| `charts/clinical.py` | Clinical slopegraphs, Shannon scatter, taxa × clinical Spearman heatmap |
| `charts/stats.py` | Wilcoxon / Mann-Whitney tests, LME trajectory, PERMANOVA (numpy-vectorised), diversity summary table |
| `charts/orchestrator.py` | `compute_chart_data()` + `ReportConfig` — coordinates all chart builders, caches Bray-Curtis matrix |
| `charts/preprocessing.py` | Shared row helpers: `get_patient_timepoints`, `sorted_timepoints`, `get_base_groups`, `filter_rows` |
| `charts/metrics.py` | `METRIC_LABELS` dict, `metric_value()`, `pielou_evenness()` — extracted to break circular import |
| `charts/stats_helpers.py` | Pure stat functions: Wilcoxon, Mann-Whitney, Welch t, Pearson r, Spearman ρ, BH-FDR |
| `charts/insights.py` | Dynamic text summaries generated from chart payloads (extracted from renderer) |
| `charts/renderer.py` | Fills HTML templates — cohort report (`template.html`) and per-patient report (`patient_template.html`) |
| `charts/template.html` | HTML/CSS/JS report shell — all controls, Plotly rendering, smooth scroll |
| `charts/patient_template.html` | Per-patient HTML report shell — stability score, radar, composition, clinical charts |
| `charts/plotly.min.js` | Bundled Plotly.js v2.35.2 — required for offline/HPC use |
| `Dockerfile` | Container image (`python:3.11-slim` + deps, no source code baked in) |
