# Changelog

All notable changes to MicroSee are documented here.  
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).  
Version numbers follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

## [0.1.0] — 2026-05-21

First public release of MicroSee. Produced for the AGB2026 course (Group D).

### Added

**Core pipeline**
- `generate_report.py` CLI entry point (`microsee-report` command via `pip install -e .`)
- QIIME2 TSV parsers: feature-table, taxonomy, metadata, alpha-diversity, distance-matrix
- Integrator: family-level aggregation, relative abundance, metadata join, alpha diversity join
- Pydantic v2 data models with full type annotations (`SampleRow`, `IntegrateResult`, etc.)

**Chart engine — 30+ interactive Plotly charts across 7 sections**
- *Taxonomy*: stacked bar (27 filter variants), top-taxa ranking, donut, sunburst
- *Alpha diversity*: strip/box/violin, significance brackets (Wilcoxon + Mann-Whitney), rarefaction curves, multi-metric panel
- *Beta diversity*: PCoA (Bray-Curtis + Jaccard), NMDS, hierarchical dendrogram, Δ abundance heatmap
- *Individual*: paired slopegraph, stability bar, diversity rank, patient radar, NMDS trajectories, small multiples
- *Comparative*: log fold-change bar, volcano (BH-FDR), ANCOM-style CLR, abundance heatmap, taxon correlation matrix
- *Clinical*: 6MWT + IL-18 slopegraphs, Shannon vs clinical scatter (Pearson r), taxa × clinical Spearman heatmap
- *Statistics*: PERMANOVA table, diversity summary table, LME-style trajectory with 95% CI

**Output modes**
- `--mode cohort` — one combined offline HTML report (~5 MB)
- `--mode patient` — one HTML per patient
- `--mode all` — cohort report + per-patient reports with cross-navigation

**Scientific transparency**
- Dynamic section insights (auto-generated text banners)
- Per-chart ℹ explanation panels (what/finding/method)
- Explicit warnings when alpha diversity is estimated from family-level data
- Explicit warnings for zero-abundance samples and unparseable timepoints
- Warnings for duplicate sample IDs (detected before pandas silently renames columns)

**Reproducibility**
- Provenance metadata embedded in every report footer: timestamp, version, git commit, Python version, platform, dependency versions, input file SHA-256 hashes
- `--provenance-json` flag to write a machine-readable JSON sidecar
- Deterministic PERMANOVA (seed derived from sorted sample IDs via MD5)
- Stable taxa ordering (sorted by mean abundance across all samples)

**Deployment**
- Self-contained offline HTML — Plotly.js v2.35.2 bundled; no CDN calls, no internet required
- HPC compatible: pure CLI, no GUI, no R, no scipy
- Nextflow process (`main.nf`) with conda, Docker, and Singularity profiles
- Dockerfile (`python:3.11-slim`, non-root `USER nobody`)
- Conda environment (`environment.yml`) with pinned deps
- PEP 561 typed package (`py.typed`)

**Testing**
- 140 tests across 7 test files
- Unit tests: parsers, chart builders, distance metrics, statistical helpers, preprocessing
- Integration tests: full HTML generation and CLI smoke tests (marked `integration`)
- End-to-end tests: reproducibility, scientific invariants, edge cases, validation
- Fixture dataset: 12 patients × 2 timepoints = 24 samples with clinical columns

**CI/CD**
- GitHub Actions: ruff lint, mypy typecheck, pytest (3.11 + 3.12), integration tests, Nextflow smoke
- Coverage report uploaded as CI artifact (Python 3.11 matrix leg)

### Fixed
- Shannon=0 false fallback: legitimately-zero diversity values were incorrectly re-estimated
- Duplicate sample IDs: raw header checked before pandas reads (pandas silently renames duplicates)
- Zero-abundance samples: now emit a warning instead of silently producing uniform abundances
- `_base_groups()` duplication: removed local copy in `comparative.py` that used `r["group"]` (KeyError risk); unified on `preprocessing.get_base_groups()`
- Dead model classes (`TSVUpload`, `ParseError`) removed from `models.py`
- UTF-8 BOM handling: `utf-8-sig` encoding in all file reads (transparent for plain UTF-8)
- Output directory creation: `mkdir(parents=True, exist_ok=True)` before writing output files

---

[Unreleased]: https://github.com/egenomics/agb2026/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/egenomics/agb2026/releases/tag/v0.1.0
