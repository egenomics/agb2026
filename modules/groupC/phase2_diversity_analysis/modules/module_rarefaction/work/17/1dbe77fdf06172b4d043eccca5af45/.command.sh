#!/bin/bash -ue
mkdir -p rarefaction_output \
             merged_tsvs \
             final_diversity/core_metrics \
             final_diversity/diversity_table/shannon \
             final_diversity/diversity_table/observed \
             final_diversity/diversity_table/faith \
             final_diversity/diversity_table/simpson \
             final_diversity/diversity_table/weighted_unifrac \
             final_diversity/diversity_table/bray_curtis

    # Read the sampling depth produced by module_diversity_metrics
    SAMPLING_DEPTH=$(cat sampling_depth.txt | tr -d '[:space:]')
    echo "Using sampling depth from module_diversity_metrics: $SAMPLING_DEPTH"

    # max_depth for the rarefaction sweep: use ext.max_depth if set, otherwise fall back to sampling depth
    MAX_DEPTH=$SAMPLING_DEPTH


    # Step 1: copy / merge tables into the working directory

    if [ "1" -gt 1 ]; then
        MERGE_ARGS=""
        for t in rarefied_table.qza; do
            MERGE_ARGS="$MERGE_ARGS --i-tables $t"
        done
        qiime feature-table merge \
            $MERGE_ARGS \
            --o-merged-table rarefaction_output/merged_table.qza
    else
        cp rarefied_table.qza rarefaction_output/merged_table.qza
    fi


    # Step 2: merge the alpha diversity TSVs from module_diversity_metrics
    # Each metric lives in its own subdirectory of diversity_table_dir,
    # so there is no filename collision when Nextflow stages the directory.

    python -c "
import pandas as pd, sys, os

TAB = chr(9)
base = 'diversity_table'

for mname, subdir in [
    ('shannon',           'shannon'),
    ('observed_features', 'observed'),
    ('faith_pd',          'faith'),
    ('simpson',           'simpson'),
]:
    fpath = os.path.join(base, subdir, 'alpha-diversity.tsv')
    try:
        df = pd.read_csv(fpath, sep=TAB, index_col=0)
        df.to_csv('merged_tsvs/' + mname + '.tsv', sep=TAB)
        print('Read ' + str(len(df)) + ' samples for ' + mname)
    except Exception as e:
        print('Warning: could not read ' + fpath + ': ' + str(e), file=sys.stderr)
"


    # Step 3: QIIME2 alpha-rarefaction curves on the merged table

    qiime diversity alpha-rarefaction \
        --i-table         rarefaction_output/merged_table.qza \
        --i-phylogeny     rooted-tree.qza \
        --p-max-depth     $MAX_DEPTH \
        --p-steps         20 \
        --p-metrics       observed_features \
        --p-metrics       shannon \
        --p-metrics       faith_pd \
        --m-metadata-file sample-metadata.tsv \
        --o-visualization rarefaction_output/rarefaction_curves.qzv


    # Step 4: Python plateau detection

    python /home/nuria/Documents/AGB/project/agb2026/modules/groupC/phase2_diversity_analysis/module_rarefaction/select_rarefaction_depth.py \
        --curves          rarefaction_output/rarefaction_curves.qzv \
        --shannon         merged_tsvs/shannon.tsv \
        --observed        merged_tsvs/observed_features.tsv \
        --faith           merged_tsvs/faith_pd.tsv \
        --simpson         merged_tsvs/simpson.tsv \
        --output          rarefaction_output \
        --method          coverage \
        --metric          observed_features \
        --coverage-pct    90 \
        --percentile      50 \
        --dropout-max     0.10 \
        --sampling-depth  $SAMPLING_DEPTH \
        

    DEPTH=$(cat rarefaction_output/rarefaction_threshold.txt | tr -d '[:space:]')
    echo "$DEPTH" > rarefaction_depth_used.txt
    echo "Batch rarefaction: auto-selected rarefaction depth = $DEPTH reads"


    # Step 5: filter samples that did not reach the plateau

    qiime feature-table filter-samples \
        --i-table          rarefaction_output/merged_table.qza \
        --m-metadata-file  rarefaction_output/samples_to_keep.tsv \
        --o-filtered-table rarefaction_output/filtered_table.qza


    # Step 6: rarefy the filtered table at the auto-selected depth

    qiime feature-table rarefy \
        --i-table          rarefaction_output/filtered_table.qza \
        --p-sampling-depth $DEPTH \
        --o-rarefied-table rarefaction_output/rarefied_table.qza


    # Step 7: core-metrics-phylogenetic on the rarefied table

    qiime diversity core-metrics-phylogenetic \
        --i-phylogeny rooted-tree.qza \
        --i-table     rarefaction_output/rarefied_table.qza \
        --p-sampling-depth $DEPTH \
        --m-metadata-file  sample-metadata.tsv \
        --o-rarefied-table                     final_diversity/core_metrics/rarefied_table.qza \
        --o-faith-pd-vector                    final_diversity/core_metrics/faith_pd_vector.qza \
        --o-observed-features-vector           final_diversity/core_metrics/observed_features_vector.qza \
        --o-shannon-vector                     final_diversity/core_metrics/shannon_vector.qza \
        --o-evenness-vector                    final_diversity/core_metrics/evenness_vector.qza \
        --o-unweighted-unifrac-distance-matrix final_diversity/core_metrics/unweighted_unifrac_distance_matrix.qza \
        --o-weighted-unifrac-distance-matrix   final_diversity/core_metrics/weighted_unifrac_distance_matrix.qza \
        --o-jaccard-distance-matrix            final_diversity/core_metrics/jaccard_distance_matrix.qza \
        --o-bray-curtis-distance-matrix        final_diversity/core_metrics/bray_curtis_distance_matrix.qza \
        --o-unweighted-unifrac-pcoa-results    final_diversity/core_metrics/unweighted_unifrac_pcoa_results.qza \
        --o-weighted-unifrac-pcoa-results      final_diversity/core_metrics/weighted_unifrac_pcoa_results.qza \
        --o-jaccard-pcoa-results               final_diversity/core_metrics/jaccard_pcoa_results.qza \
        --o-bray-curtis-pcoa-results           final_diversity/core_metrics/bray_curtis_pcoa_results.qza \
        --o-unweighted-unifrac-emperor         final_diversity/core_metrics/unweighted_unifrac_emperor.qzv \
        --o-weighted-unifrac-emperor           final_diversity/core_metrics/weighted_unifrac_emperor.qzv \
        --o-jaccard-emperor                    final_diversity/core_metrics/jaccard_emperor.qzv \
        --o-bray-curtis-emperor                final_diversity/core_metrics/bray_curtis_emperor.qzv

    qiime diversity alpha \
        --i-table   rarefaction_output/rarefied_table.qza \
        --p-metric  simpson \
        --o-alpha-diversity final_diversity/core_metrics/simpson_vector.qza


    # Step 8: export final diversity metrics to TSV

    qiime tools export \
        --input-path  final_diversity/core_metrics/shannon_vector.qza \
        --output-path final_diversity/diversity_table/shannon

    qiime tools export \
        --input-path  final_diversity/core_metrics/observed_features_vector.qza \
        --output-path final_diversity/diversity_table/observed

    qiime tools export \
        --input-path  final_diversity/core_metrics/faith_pd_vector.qza \
        --output-path final_diversity/diversity_table/faith

    qiime tools export \
        --input-path  final_diversity/core_metrics/simpson_vector.qza \
        --output-path final_diversity/diversity_table/simpson

    qiime tools export \
        --input-path  final_diversity/core_metrics/weighted_unifrac_distance_matrix.qza \
        --output-path final_diversity/diversity_table/weighted_unifrac

    qiime tools export \
        --input-path  final_diversity/core_metrics/bray_curtis_distance_matrix.qza \
        --output-path final_diversity/diversity_table/bray_curtis
