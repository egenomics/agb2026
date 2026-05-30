// Module 2: RAREFACTION_THRESHOLD

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
          path(diversity_table_dir),   // entire diversity_table/ folder: avoids alpha-diversity.tsv name collision
          path(sampling_depth_file)

    output:
    tuple val(meta), path("rarefaction_output/rarefaction_threshold.txt"),       emit: threshold
    tuple val(meta), path("rarefaction_output/rarefied_table.qza"),              emit: rarefied_table
    tuple val(meta), path("rarefaction_output/rarefaction_curves.qzv"),          emit: curves
    tuple val(meta), path("rarefaction_output/sample_qc.tsv"),                   emit: qc
    tuple val(meta), path("rarefaction_output/report.txt"),                      emit: report
    tuple val(meta), path("rarefaction_output/rarefaction_plots.png"),              emit: plots
    tuple val(meta), path("final_diversity/diversity_table/shannon/*"),          emit: shannon
    tuple val(meta), path("final_diversity/diversity_table/observed/*"),         emit: observed_features
    tuple val(meta), path("final_diversity/diversity_table/faith/*"),            emit: faith_pd
    tuple val(meta), path("final_diversity/diversity_table/simpson/*"),          emit: simpson
    tuple val(meta), path("final_diversity/diversity_table/weighted_unifrac/*"), emit: weighted_unifrac
    tuple val(meta), path("final_diversity/diversity_table/bray_curtis/*"),      emit: bray_curtis
    tuple val(meta), path("final_diversity/core_metrics/*.qza"),                 emit: qza_files
    tuple val(meta), path("final_diversity/core_metrics/*.qzv"),                 emit: qzv_files
    tuple val(meta), path("rarefaction_depth_used.txt"),                         emit: depth_used

    when:
    task.ext.when == null || task.ext.when

    script:
    def method       = task.ext.method       ?: 'coverage'
    def div_metric   = task.ext.metric       ?: 'observed_features'
    def coverage_pct = task.ext.coverage_pct ?: 90
    def percentile   = task.ext.percentile   ?: 50
    def dropout_max  = task.ext.dropout_max  ?: 0.10
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

    # Read the sampling depth produced by module_diversity_metrics
    SAMPLING_DEPTH=\$(cat ${sampling_depth_file} | tr -d '[:space:]')

    # max_depth for the rarefaction sweep: use ext.max_depth if set, otherwise fall back to sampling depth
    MAX_DEPTH=${task.ext.max_depth ?: '\$SAMPLING_DEPTH'}


    # Step 1: copy / merge tables into the working directory

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


    # Step 2: merge the alpha diversity TSVs from module_diversity_metrics
    # Each metric lives in its own subdirectory of diversity_table dir,
    # so there is no filename collision when Nextflow stages the directory.

    python -c "
import pandas as pd, sys, os

TAB = chr(9)
base = '${diversity_table_dir}'

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

    qiime diversity alpha-rarefaction \\
        --i-table         rarefaction_output/merged_table.qza \\
        --i-phylogeny     ${phylogeny} \\
        --p-max-depth     \$MAX_DEPTH \\
        --p-steps         ${steps} \\
        --p-metrics       observed_features \\
        --p-metrics       shannon \\
        --p-metrics       faith_pd \\
        --m-metadata-file ${metadata} \\
        --o-visualization rarefaction_output/rarefaction_curves.qzv


    # Step 4: Python plateau detection

    python ${projectDir}/select_rarefaction_depth.py \\
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
        --sampling-depth  \$SAMPLING_DEPTH \\
        ${args}

    DEPTH=\$(cat rarefaction_output/rarefaction_threshold.txt | tr -d '[:space:]')
    echo "\$DEPTH" > rarefaction_depth_used.txt

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

    """
}

process SKIP_RAREFACTION {
    tag "$meta.id"

    conda params.qiime2_conda_env ?: "qiime2=2026.1"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qiime2-amplicon-2024.10:latest' :
        'quay.io/qiime2/amplicon:2024.10' }"

    input:
    tuple val(meta), val(n_samples), val(min_batch_size)

    output:
    tuple val(meta), path("rarefaction_skipped_report.png"), emit: report

    script:
    """
    python << 'PYTHON'
from PIL import Image, ImageDraw, ImageFont
from datetime import datetime

width = 300
height = 200

img = Image.new("RGB", (width, height), "white")
draw = ImageDraw.Draw(img)

try:
    font_title = ImageFont.truetype("DejaVuSansMono.ttf", 32)
    font_text  = ImageFont.truetype("DejaVuSansMono.ttf", 24)
except:
    font_title = ImageFont.load_default()
    font_text  = ImageFont.load_default()

lines = [
    "RAREFACTION SKIPPED",
    "",
    f"Sample count    : ${n_samples}",
    f"Minimum required: ${min_batch_size}",
    "",
    "Automatic rarefaction threshold selection was not",
    f"performed because the number of samples (${n_samples})",
    f"did not reach the minimum required (${min_batch_size}).",
    "",
    "Action required: accumulate more samples and",
    "re-run the pipeline once the threshold is met."
]

y = 20

draw.text((20, y), lines[0], fill="black", font=font_title)
y += 20

for line in lines[2:]:
    draw.text((20, y), line, fill="black", font=font_text)
    y += 15

img.save("rarefaction_skipped_report.png")
PYTHON
    """
}


// Workflow: 

workflow {

    def metrics_dir = params.diversity_metrics_dir

    def rarefied_table      = file("${metrics_dir}/core_metrics_output/rarefied_table.qza")
    def phylogeny           = file("${metrics_dir}/data/rooted-tree.qza")
    def metadata            = file("${metrics_dir}/data/sample-metadata.tsv")
    def sampling_depth_f    = file("${metrics_dir}/sampling_depth.txt")
    def diversity_table_dir = file("${metrics_dir}/diversity_table")

    // Count samples in each alpha diversity TSV and take the minimum.
    // The minimum is what rarefaction will actually operate on — if one metric
    // has fewer samples (e.g. simpson calculated on unrarefied table), the min
    // correctly reflects the most restrictive case.
    def tsv_subpaths = ['shannon/alpha-diversity.tsv',
                        'observed/alpha-diversity.tsv',
                        'faith/alpha-diversity.tsv',
                        'simpson/alpha-diversity.tsv']

    def per_metric_counts = tsv_subpaths.collect { subpath ->
        def tsv = file("${diversity_table_dir}/${subpath}")
        if (tsv.exists()) {
            tsv.readLines().findAll { !it.startsWith('#') && !it.trim().isEmpty() }.size() - 1
        } else {
            null
        }
    }.findAll { it != null }

    def n_samples = per_metric_counts.min()

    log.info "Sample counts per metric: ${tsv_subpaths.collect { it.split('/')[0] }.join(', ')} = ${per_metric_counts.join(', ')} → min = ${n_samples}"

    def min_size = params.min_batch_size
    def meta     = [id: "batch_n${n_samples}"]

    // Branch: run rarefaction or emit a skip report
    if (n_samples >= min_size) {

        def input_ch = Channel.of(
            tuple(
                meta,
                rarefied_table,
                phylogeny,
                metadata,
                diversity_table_dir,
                sampling_depth_f
            )
        )
        RAREFACTION_THRESHOLD(input_ch)

    } else {

        log.warn "Only ${n_samples} samples found — minimum is ${min_size}. Generating skip report."

        SKIP_RAREFACTION(
            Channel.of(tuple(meta, n_samples, min_size))
        )

        SKIP_RAREFACTION.out.report.view { meta_val, report ->
            "Rarefaction skipped — report written to: ${report}"
        }
    }
}
