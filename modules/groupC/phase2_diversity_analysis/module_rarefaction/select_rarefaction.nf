// =============================================================================
// Module 2: RAREFACTION_THRESHOLD
// =============================================================================

process RAREFACTION_THRESHOLD {
    tag "$meta.id"
    label 'process_medium'

    conda params.qiime2_conda_env ?: "qiime2=2026.1"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qiime2-amplicon-2024.10:latest' :
        'quay.io/qiime2/amplicon:2024.10' }"

    input:
    tuple val(meta),
          path(tables),
          path(phylogeny),
          path(metadata),
          path(shannon_tsvs),
          path(observed_tsvs),
          path(faith_tsvs),
          path(simpson_tsvs)

    output:
    tuple val(meta), path("rarefaction_output/rarefaction_threshold.txt"),       emit: threshold
    tuple val(meta), path("rarefaction_output/rarefied_table.qza"),              emit: rarefied_table
    tuple val(meta), path("rarefaction_output/rarefaction_curves.qzv"),          emit: curves
    tuple val(meta), path("rarefaction_output/sample_qc.tsv"),                   emit: qc
    tuple val(meta), path("rarefaction_output/report.txt"),                      emit: report
    tuple val(meta), path("rarefaction_output/*.pdf"),                           emit: plots
    tuple val(meta), path("final_diversity/diversity_table/shannon/*"),          emit: shannon
    tuple val(meta), path("final_diversity/diversity_table/observed/*"),         emit: observed_features
    tuple val(meta), path("final_diversity/diversity_table/faith/*"),            emit: faith_pd
    tuple val(meta), path("final_diversity/diversity_table/simpson/*"),          emit: simpson
    tuple val(meta), path("final_diversity/diversity_table/weighted_unifrac/*"), emit: weighted_unifrac
    tuple val(meta), path("final_diversity/diversity_table/bray_curtis/*"),      emit: bray_curtis
    tuple val(meta), path("final_diversity/core_metrics/*.qza"),                 emit: qza_files
    tuple val(meta), path("final_diversity/core_metrics/*.qzv"),                 emit: qzv_files
    tuple val(meta), path("rarefaction_depth_used.txt"),                         emit: depth_used
    path "versions.yml",                                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def method       = task.ext.method       ?: 'coverage'
    def div_metric   = task.ext.metric       ?: 'observed_features'  

    def coverage_pct = task.ext.coverage_pct ?: 90
    def percentile   = task.ext.percentile   ?: 50
    def dropout_max  = task.ext.dropout_max  ?: 0.10
    def max_depth = task.ext.max_depth ?: params.sampling_depth
                                                                    
    def steps        = task.ext.steps        ?: 20
    def args         = task.ext.args         ?: ''

    def table_count  = tables instanceof List ? tables.size() : 1

    """
    mkdir -p rarefaction_output \\
             merged_tsvs \\
             final_diversity/core_metrics \\
             final_diversity/diversity_table/shannon \\
             final_diversity/diversity_table/observed \\
             final_diversity/diversity_table/faith \\
             final_diversity/diversity_table/simpson \\
             final_diversity/diversity_table/weighted_unifrac \\
             final_diversity/diversity_table/bray_curtis


    # Step 1: merge all tables in the batch into one

    if [ "${table_count}" -gt 1 ]; then
        MERGE_ARGS=""
        for t in ${tables}; do
            MERGE_ARGS="\$MERGE_ARGS --i-tables \$t"
        done
        qiime feature-table merge \\
            \$MERGE_ARGS \\
            --o-merged-table rarefaction_output/merged_table.qza
    else
        cp ${tables} rarefaction_output/merged_table.qza
    fi


    # Step 2: merge the alpha diversity TSVs from Module 1

    python -c "
import pandas as pd, sys

TAB = chr(9)

for mname, flist in [
    ('shannon', '${shannon_tsvs}'.split()),
    ('observed_features', '${observed_tsvs}'.split()),
    ('faith_pd', '${faith_tsvs}'.split()),
    ('simpson', '${simpson_tsvs}'.split()),
]:
    frames = []
    for fpath in flist:
        try:
            frames.append(pd.read_csv(fpath, sep=TAB, index_col=0))
        except Exception as e:
            print('Warning: could not read ' + fpath + ': ' + str(e), file=sys.stderr)
    if frames:
        combined = pd.concat(frames)
        combined.to_csv('merged_tsvs/' + mname + '.tsv', sep=TAB)
        print('Merged ' + str(len(frames)) + ' files for ' + mname + ': ' + str(len(combined)) + ' samples')
"


    # Step 3: QIIME2 alpha-rarefaction curves on the merged table

    qiime diversity alpha-rarefaction \\
        --i-table         rarefaction_output/merged_table.qza \\
        --i-phylogeny     ${phylogeny} \\
        --p-max-depth     ${max_depth} \\
        --p-steps         ${steps} \\
        --p-metrics       observed_features \\
        --p-metrics       shannon \\
        --p-metrics       faith_pd \\
        --m-metadata-file ${metadata} \\
        --o-visualization rarefaction_output/rarefaction_curves.qzv


    # Step 4: Python plateau detection 


    python ${projectDir}/bin/select_rarefaction_depth.py \\
        --curves          rarefaction_output/rarefaction_curves.qzv \\
        --shannon         merged_tsvs/shannon.tsv \\
        --observed        merged_tsvs/observed_features.tsv \\
        --faith           merged_tsvs/faith_pd.tsv \\
        --simpson         merged_tsvs/simpson.tsv \\
        --output          rarefaction_output \\
        --method          ${method} \\
        --metric          ${div_metric} \\
        --coverage-pct    ${coverage_pct} \\
        --percentile      ${percentile} \\
        --dropout-max     ${dropout_max} \\
        --sampling-depth  ${params.sampling_depth} \\
        ${args}

    DEPTH=\$(cat rarefaction_output/rarefaction_threshold.txt | tr -d '[:space:]')
    echo "\$DEPTH" > rarefaction_depth_used.txt
    echo "Batch ${meta.id}: auto-selected rarefaction depth = \$DEPTH reads"


    # Step 5: filter samples that did not reach the plateau

    qiime feature-table filter-samples \\
        --i-table          rarefaction_output/merged_table.qza \\
        --m-metadata-file  rarefaction_output/samples_to_keep.tsv \\
        --o-filtered-table rarefaction_output/filtered_table.qza


    # Step 6: rarefy the filtered table at the auto-selected depth

    qiime feature-table rarefy \\
        --i-table          rarefaction_output/filtered_table.qza \\
        --p-sampling-depth \$DEPTH \\
        --o-rarefied-table rarefaction_output/rarefied_table.qza


    # Step 7: core-metrics-phylogenetic on the rarefied table
    # 
    qiime diversity core-metrics-phylogenetic \\
        --i-phylogeny ${phylogeny} \\
        --i-table     rarefaction_output/rarefied_table.qza \\
        --p-sampling-depth \$DEPTH \\
        --m-metadata-file  ${metadata} \\
        --o-rarefied-table                     final_diversity/core_metrics/rarefied_table.qza \\
        --o-faith-pd-vector                    final_diversity/core_metrics/faith_pd_vector.qza \\
        --o-observed-features-vector           final_diversity/core_metrics/observed_features_vector.qza \\
        --o-shannon-vector                     final_diversity/core_metrics/shannon_vector.qza \\
        --o-evenness-vector                    final_diversity/core_metrics/evenness_vector.qza \\
        --o-unweighted-unifrac-distance-matrix final_diversity/core_metrics/unweighted_unifrac_distance_matrix.qza \\
        --o-weighted-unifrac-distance-matrix   final_diversity/core_metrics/weighted_unifrac_distance_matrix.qza \\
        --o-jaccard-distance-matrix            final_diversity/core_metrics/jaccard_distance_matrix.qza \\
        --o-bray-curtis-distance-matrix        final_diversity/core_metrics/bray_curtis_distance_matrix.qza \\
        --o-unweighted-unifrac-pcoa-results    final_diversity/core_metrics/unweighted_unifrac_pcoa_results.qza \\
        --o-weighted-unifrac-pcoa-results      final_diversity/core_metrics/weighted_unifrac_pcoa_results.qza \\
        --o-jaccard-pcoa-results               final_diversity/core_metrics/jaccard_pcoa_results.qza \\
        --o-bray-curtis-pcoa-results           final_diversity/core_metrics/bray_curtis_pcoa_results.qza \\
        --o-unweighted-unifrac-emperor         final_diversity/core_metrics/unweighted_unifrac_emperor.qzv \\
        --o-weighted-unifrac-emperor           final_diversity/core_metrics/weighted_unifrac_emperor.qzv \\
        --o-jaccard-emperor                    final_diversity/core_metrics/jaccard_emperor.qzv \\
        --o-bray-curtis-emperor                final_diversity/core_metrics/bray_curtis_emperor.qzv

    qiime diversity alpha \\
        --i-table   rarefaction_output/rarefied_table.qza \\
        --p-metric  simpson \\
        --o-alpha-diversity final_diversity/core_metrics/simpson_vector.qza


    # Step 8: export final diversity metrics to TSV

    qiime tools export \\
        --input-path  final_diversity/core_metrics/shannon_vector.qza \\
        --output-path final_diversity/diversity_table/shannon

    qiime tools export \\
        --input-path  final_diversity/core_metrics/observed_features_vector.qza \\
        --output-path final_diversity/diversity_table/observed

    qiime tools export \\
        --input-path  final_diversity/core_metrics/faith_pd_vector.qza \\
        --output-path final_diversity/diversity_table/faith

    qiime tools export \\
        --input-path  final_diversity/core_metrics/simpson_vector.qza \\
        --output-path final_diversity/diversity_table/simpson

    qiime tools export \\
        --input-path  final_diversity/core_metrics/weighted_unifrac_distance_matrix.qza \\
        --output-path final_diversity/diversity_table/weighted_unifrac

    qiime tools export \\
        --input-path  final_diversity/core_metrics/bray_curtis_distance_matrix.qza \\
        --output-path final_diversity/diversity_table/bray_curtis

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$(qiime info 2>&1 | grep "QIIME 2 release" | sed 's/.*: //')
        python: \$(python --version 2>&1 | sed 's/Python //')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
