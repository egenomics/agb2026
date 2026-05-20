#!/bin/bash
# test script to run the rarefaction module manually (without nextflow)

set -euo pipefail # stop if anything fails


# paths
TABLES="/data/AGB/project/agb2026/modules/module_diversity_metrics/results/diversity_table/core_metrics_output/rarefied_table.qza"
PHYLOGENY="/data/AGB/project/agb2026/modules/module_diversity_metrics/data/rooted-tree.qza"
METADATA="/data/AGB/project/agb2026/modules/module_diversity_metrics/data/sample-metadata.tsv"

# alpha diversity outputs from module 1
SHANNON_TSVS="/data/AGB/project/agb2026/modules/module_diversity_metrics/results/diversity_table/diversity_table/shannon/alpha-diversity.tsv"
OBSERVED_TSVS="/data/AGB/project/agb2026/modules/module_diversity_metrics/results/diversity_table/diversity_table/observed/alpha-diversity.tsv"
FAITH_TSVS="/data/AGB/project/agb2026/modules/module_diversity_metrics/results/diversity_table/diversity_table/faith/alpha-diversity.tsv"
SIMPSON_TSVS="/data/AGB/project/agb2026/modules/module_diversity_metrics/results/diversity_table/diversity_table/simpson/alpha-diversity.tsv"


# parameters

# which method to use for picking the rarefaction depth
METHOD="coverage"

# for coverage method: percentage of samples that should have plateaued
COVERAGE_PCT="90"

# which metric to use for plateau detection
METRIC="observed_features"

# the depth that was used in module 1
SAMPLING_DEPTH="1103"

# max depth for the rarefaction curve sweep
MAX_DEPTH="1000"

STEPS="20"
DROPOUT_MAX="0.10"

PROJECT_DIR="/data/AGB/project/agb2026"


# make output folders

mkdir -p rarefaction_output \
         merged_tsvs \
         final_diversity/core_metrics \
         final_diversity/diversity_table/shannon \
         final_diversity/diversity_table/observed \
         final_diversity/diversity_table/faith \
         final_diversity/diversity_table/simpson \
         final_diversity/diversity_table/weighted_unifrac \
         final_diversity/diversity_table/bray_curtis


# step 1: copy the table into the working directory
 
cp "$TABLES" rarefaction_output/merged_table.qza

 # step 2: merge the alpha diversity tsvs from module 1
 
# combine the files
python -c "
import pandas as pd, sys

for metric, pattern in [
    ('shannon', '$SHANNON_TSVS'.split()),
    ('observed_features', '$OBSERVED_TSVS'.split()),
    ('faith_pd', '$FAITH_TSVS'.split()),
    ('simpson', '$SIMPSON_TSVS'.split()),
]:
    frames = []
    for f in pattern:
        try:
            frames.append(pd.read_csv(f, sep='\t', index_col=0))
        except Exception as e:
            print(f'Warning: could not read {f}: {e}', file=sys.stderr)
    if frames:
        combined = pd.concat(frames)
        combined.to_csv(f'merged_tsvs/{metric}.tsv', sep='\t')
        print(f'Merged {len(frames)} files for {metric}: {len(combined)} samples')
"

 
# step 3: run qiime2 to generate the rarefaction curves
 
qiime diversity alpha-rarefaction \
    --i-table rarefaction_output/merged_table.qza \
    --i-phylogeny "$PHYLOGENY" \
    --p-max-depth "$MAX_DEPTH" \
    --p-steps "$STEPS" \
    --p-metrics observed_features \
    --p-metrics shannon \
    --p-metrics faith_pd \
    --m-metadata-file "$METADATA" \
    --o-visualization rarefaction_output/rarefaction_curves.qzv

 
# step 4: run the plateau detection script to pick a threshold
 
python "$PROJECT_DIR/modules/module_rarefaction/select_rarefaction_depth.py" \
    --curves rarefaction_output/rarefaction_curves.qzv \
    --shannon merged_tsvs/shannon.tsv \
    --observed merged_tsvs/observed_features.tsv \
    --faith merged_tsvs/faith_pd.tsv \
    --simpson merged_tsvs/simpson.tsv \
    --output rarefaction_output \
    --method "$METHOD" \
    --metric "$METRIC" \
    --coverage-pct "$COVERAGE_PCT" \
    --dropout-max "$DROPOUT_MAX" \
    --sampling-depth "$SAMPLING_DEPTH"

# read the chosen depth from the output file
DEPTH=$(cat rarefaction_output/rarefaction_threshold.txt | tr -d '[:space:]')
echo "$DEPTH" > rarefaction_depth_used.txt
echo "selected depth: $DEPTH reads"

 
# step 5: drop samples that are below the threshold
 
qiime feature-table filter-samples \
    --i-table rarefaction_output/merged_table.qza \
    --m-metadata-file rarefaction_output/samples_to_keep.tsv \
    --o-filtered-table rarefaction_output/filtered_table.qza

 
# step 6: rarefy to the selected depth
 
qiime feature-table rarefy \
    --i-table rarefaction_output/filtered_table.qza \
    --p-sampling-depth "$DEPTH" \
    --o-rarefied-table rarefaction_output/rarefied_table.qza

echo "done! rarefied table is at rarefaction_output/rarefied_table.qza"
