#!/bin/bash -ue
mkdir -p core_metrics_output \
  diversity_table/shannon \
         diversity_table/observed \
         diversity_table/faith \
         diversity_table/simpson \
         diversity_table/weighted_unifrac \
         diversity_table/bray_curtis

qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree.qza \
    --i-table table.qza \
    --p-sampling-depth 1103 \
    --m-metadata-file sample-metadata.tsv \
    --o-rarefied-table                     core_metrics_output/rarefied_table.qza \
    --o-faith-pd-vector                    core_metrics_output/faith_pd_vector.qza \
    --o-observed-features-vector           core_metrics_output/observed_features_vector.qza \
    --o-shannon-vector                     core_metrics_output/shannon_vector.qza \
    --o-evenness-vector                    core_metrics_output/evenness_vector.qza \
    --o-unweighted-unifrac-distance-matrix core_metrics_output/unweighted_unifrac_distance_matrix.qza \
    --o-weighted-unifrac-distance-matrix   core_metrics_output/weighted_unifrac_distance_matrix.qza \
    --o-jaccard-distance-matrix            core_metrics_output/jaccard_distance_matrix.qza \
    --o-bray-curtis-distance-matrix        core_metrics_output/bray_curtis_distance_matrix.qza \
    --o-unweighted-unifrac-pcoa-results    core_metrics_output/unweighted_unifrac_pcoa_results.qza \
    --o-weighted-unifrac-pcoa-results      core_metrics_output/weighted_unifrac_pcoa_results.qza \
    --o-jaccard-pcoa-results               core_metrics_output/jaccard_pcoa_results.qza \
    --o-bray-curtis-pcoa-results           core_metrics_output/bray_curtis_pcoa_results.qza \
    --o-unweighted-unifrac-emperor         core_metrics_output/unweighted_unifrac_emperor.qzv \
    --o-weighted-unifrac-emperor           core_metrics_output/weighted_unifrac_emperor.qzv \
    --o-jaccard-emperor                    core_metrics_output/jaccard_emperor.qzv \
    --o-bray-curtis-emperor                core_metrics_output/bray_curtis_emperor.qzv

qiime diversity alpha \
    --i-table table.qza \
    --p-metric simpson \
    --o-alpha-diversity core_metrics_output/simpson_vector.qza

qiime tools export \
    --input-path core_metrics_output/shannon_vector.qza \
    --output-path diversity_table/shannon

qiime tools export \
    --input-path core_metrics_output/observed_features_vector.qza \
    --output-path diversity_table/observed

qiime tools export \
    --input-path core_metrics_output/faith_pd_vector.qza \
    --output-path diversity_table/faith

qiime tools export \
    --input-path core_metrics_output/simpson_vector.qza \
    --output-path diversity_table/simpson

qiime tools export \
    --input-path core_metrics_output/weighted_unifrac_distance_matrix.qza \
    --output-path diversity_table/weighted_unifrac

qiime tools export \
    --input-path core_metrics_output/bray_curtis_distance_matrix.qza \
    --output-path diversity_table/bray_curtis
