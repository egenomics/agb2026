#!/usr/bin/bash
set -euo pipefail

data_dir="$1"
outdir="$2"

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "$data_dir/rooted-tree.qza" \
  --i-table "$data_dir/table.qza" \
  --p-sampling-depth 1103 \
  --m-metadata-file "$data_dir/sample-metadata.tsv" \
  --output-dir "$outdir"

qiime diversity alpha \
  --i-table "$data_dir/table.qza" \
  --p-metric simpson \
  --o-alpha-diversity "$outdir/simpson_vector.qza"

qiime tools export \
  --input-path "$outdir/shannon_vector.qza" \
  --output-path diversity_table/shannon

qiime tools export \
  --input-path "$outdir/observed_features_vector.qza" \
  --output-path diversity_table/observed

qiime tools export \
  --input-path "$outdir/faith_pd_vector.qza" \
  --output-path diversity_table/faith

qiime tools export \
  --input-path "$outdir/simpson_vector.qza" \
  --output-path diversity_table/simpson

qiime tools export \
  --input-path "$outdir/weighted_unifrac_distance_matrix.qza" \
  --output-path diversity_table/weighted_unifrac



