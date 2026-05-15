# MicroSee

Interactive web application for exploratory analysis of gut microbiome data from clinical intervention studies. Upload a metadata TSV and a taxonomy TSV and get an instant, publication-ready dashboard of charts — no coding required.

---

## Features

- **Alpha diversity** — Shannon H′, Simpson, Pielou J′, Chao1, Observed Taxa, rarefaction curves
- **Beta diversity** — PCoA (Bray-Curtis & Jaccard), NMDS individual trajectories, hierarchical clustering dendrogram, Δ abundance heatmap
- **Taxonomy** — stacked bar plots, donut charts per group, taxon correlation matrix, abundance heatmap
- **Individual / longitudinal** — paired slopegraph, stability score (Bray-Curtis), diversity rank plot, patient radar profile, faceted small multiples, longitudinal diversity curves
- **Clinical** — 6-min walk test and IL-18 cytokine individual trajectories, correlation with diversity metrics
- **Statistics** — Volcano plot (Welch t-test, Benjamini-Hochberg FDR), differential abundance (log₂ fold change), PERMANOVA (99 permutations), diversity summary table

---

## Project structure

```
microsee/
├── microsee_backend/        # FastAPI — file parsing only
│   ├── app/
│   │   ├── main.py
│   │   ├── models/
│   │   ├── routes/          # POST /parse  →  SampleRow[]
│   │   └── services/        # TSV parsing & alpha-diversity calculation
│   ├── tests/
│   ├── conftest.py
│   └── requirements.txt
└── microsee_frontend/       # React + TypeScript + Vite
    ├── src/
    │   ├── charts/
    │   │   ├── alpha/       # MultiMetricAlpha, RarefactionCurves
    │   │   ├── beta/        # PCoAPlot, DendrogramPlot, DeltaHeatmap, distances.ts
    │   │   ├── clinical/    # ClinicalCharts (slopegraph + correlation)
    │   │   ├── comparative/ # ComparativeCharts (volcano, differential, heatmap, correlation)
    │   │   ├── individual/  # IndividualCharts, NMDSTrajectories, PatientCharts
    │   │   ├── longitudinal/# LongitudinalStats (chart + PERMANOVA table)
    │   │   ├── shared/      # usePlotly.ts — shared layout, colours, statistics helpers
    │   │   └── taxonomy/    # StackedBar, DonutPerGroup
    │   ├── components/      # ChartCard, FileUploader, ExpandModal
    │   ├── pages/           # ResultsPage, UploadPage
    │   ├── store/           # Zustand store
    │   ├── types/           # SampleRow type
    │   └── services/        # API client
    ├── index.html
    ├── package.json
    ├── tsconfig.json
    └── vite.config.ts
```

---

## Getting started

### Backend

```bash
cd microsee_backend
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

### Frontend

```bash
cd microsee_frontend
npm install
npm run dev        # http://localhost:5173
```

### Production build

```bash
cd microsee_frontend
npm run build      # output in dist/
```

---

## Input format

Upload two TSV files on the upload page.

### `metadata.tsv` (required)

| column | description |
|--------|-------------|
| `sample_id` | unique sample identifier |
| `group` | treatment/timepoint subgroup, e.g. `PRO_T0` |
| `time` | numeric timepoint in days (e.g. `0`, `84`) |
| `patient` | patient ID (optional — inferred from `sample_id` if absent) |
| `base_group` | baseline group label, e.g. `PRO` (optional) |
| `shannon` | Shannon H′ diversity index |
| `simpson` | Simpson diversity index |
| `chao1` | Chao1 richness estimator |
| `sixmwt` | 6-minute walk test distance in metres (optional) |
| `il18` | IL-18 cytokine level in pg/mL (optional) |
| `nmds1`, `nmds2` | pre-computed NMDS coordinates (optional) |

### `taxonomy.tsv` (required)

Rows = samples, columns = family-level taxa. Values are relative abundances (any scale — they are normalised internally).

| `sample_id` | Lachnospiraceae | Ruminococcaceae | … |
|-------------|-----------------|-----------------|---|
| S01_T0      | 28.4            | 17.2            | … |

---

## Tech stack

| Layer | Technology |
|-------|-----------|
| Frontend framework | React 18 + TypeScript |
| Build tool | Vite 5 |
| Charts | Plotly.js via react-plotly.js |
| State management | Zustand |
| Data fetching | TanStack Query |
| Backend | FastAPI + Uvicorn |
| Data processing | pandas + NumPy |

---

## Statistical methods

- **Alpha diversity**: Shannon H′, Simpson (1 − D), Pielou J′, Chao1 (F1²/2F2), Observed Taxa
- **Beta diversity**: Bray-Curtis dissimilarity, Jaccard (5% threshold), PCoA via classical MDS, UPGMA hierarchical clustering
- **Differential abundance**: Welch t-test with Satterthwaite degrees of freedom, exact p-values via regularised incomplete beta, Benjamini-Hochberg FDR, log₂ fold change with pseudo-count
- **Correlation**: Spearman rank correlation (diversity vs clinical variables)
- **PERMANOVA**: 99 permutations, Bray-Curtis distances, data-derived seed

## Mermaid Workflow

```mermaid
graph TD
    %% Data Ingestion Phase
    Input[/Input Files/] --> Parsers[Data Parsers Script]
    Parsers --> Storage[(Global Storage)]

    %% Storage Details
    subgraph StorageStructure [Storage Contents]
    ID[Sample ID]
    Meta[Metadata: Case/Control, Date, Age, Batch]
    Values[Actual Data Values]
    end
    Storage --- StorageStructure

    %% Branching Logic
    Storage --> S1[Script 1: Global Aggregator]
    Storage --> S2[Script 2: Filtered Generator]

    %% Output Phase
    S1 --> HTML1([Interactive .html Global Report])
    
    S2 --> S2Input{User Selection}
    S2Input -->|By ID| HTML2([Filtered .html Report])
    S2Input -->|By Metadata| HTML2
    ```
