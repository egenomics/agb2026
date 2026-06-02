# Group D — pipeline inputs

Map upstream QIIME2 exports to `microsee-report` / Nextflow parameters.

| CLI / Nextflow param | QIIME2 artifact | Notes |
|----------------------|-----------------|-------|
| `--feature-table` / `params.feature_table` | `feature-table.tsv` (BIOM → TSV) | ASV counts; first column = feature ID |
| `--taxonomy` / `params.taxonomy` | `taxonomy.tsv` | `Feature ID` + semicolon lineage (`f__Family`) |
| `--metadata` / `params.metadata` | `sample-metadata.tsv` | `sample-id`, group, timepoint, optional `sixmwt` / `il18` |
| `--alpha` / `params.alpha` | Merged alpha TSV | **Recommended.** Merge per-metric exports into one file |
| `--distance-matrix` / `params.distance_matrix` | `distance-matrix.tsv` | Optional Bray-Curtis (or other) matrix from QIIME2 |
| `params.mode` | — | `cohort` (default), `patient`, or `all` |
| `params.outdir` | — | Published HTML output directory |

## Alpha diversity merge (QIIME2)

QIIME2 exports one metric per file. Merge before reporting:

```bash
# Example: join on sample-id (adjust filenames to your exports)
paste -d '\t' \
  <(cut -f1,2 shannon_vector.tsv) \
  <(cut -f2 simpson_vector.tsv) \
  <(cut -f2 observed_features.tsv) \
  <(cut -f2 faith_pd_vector.tsv) \
  > alpha-diversity.tsv
```

Without `--alpha`, Shannon and Simpson are estimated from **family-level** abundances and will be lower than ASV-level QIIME2 metrics.

## Sample ID checklist

- IDs must match exactly across feature table, taxonomy (feature IDs only), metadata, alpha, and distance matrix.
- Strip trailing whitespace: `sed -i 's/[[:space:]]*$//' metadata.tsv`
- Use one separator style (`_` vs `-`) consistently.

## Nextflow examples

```bash
# Bundled fixtures (smoke test)
nextflow run workflows/groupD.nf -profile test,conda

# Production (conda — recommended until GHCR image is available)
nextflow run workflows/groupD.nf -profile conda \
  --feature_table results/feature-table.tsv \
  --taxonomy      results/taxonomy.tsv \
  --metadata      results/metadata.tsv \
  --alpha         results/alpha-diversity.tsv \
  --outdir        results/groupD/

# With precomputed distance matrix
nextflow run workflows/groupD.nf -profile conda \
  --feature_table ... \
  --distance_matrix results/bray-curtis-distance-matrix.tsv \
  ...
```
