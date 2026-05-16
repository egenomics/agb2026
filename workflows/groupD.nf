/*
 * Group D — Visualisation workflow
 *
 * Takes QIIME2 TSV exports from upstream groups and produces a single
 * self-contained HTML report (microsee_report.html) that can be opened
 * in any browser with no server required.
 *
 * Required params (set in nextflow.config or via --param on the CLI):
 *   params.feature_table  path to feature-table.tsv
 *   params.taxonomy       path to taxonomy.tsv
 *   params.metadata       path to metadata.tsv
 *   params.alpha          path to alpha-diversity.tsv  (optional, omit or set to '')
 *   params.outdir         where to publish the HTML report
 */

include { MICROSEE_REPORT } from '../modules/groupD/microsee_report/main'

workflow GROUPD {

    take:
    feature_table   // channel: path
    taxonomy        // channel: path
    metadata        // channel: path
    alpha           // channel: path  (can be file('NO_FILE'))

    main:
    MICROSEE_REPORT(
        Channel.of(file("${projectDir}/modules/groupD/microsee_report/report_generator")),
        feature_table,
        taxonomy,
        metadata,
        alpha,
    )

    emit:
    report = MICROSEE_REPORT.out.report
}

/*
 * Standalone entrypoint — run this workflow directly:
 *
 *   nextflow run workflows/groupD.nf \
 *       --feature_table path/to/feature-table.tsv \
 *       --taxonomy      path/to/taxonomy.tsv       \
 *       --metadata      path/to/metadata.tsv       \
 *       [--alpha        path/to/alpha-diversity.tsv] \
 *       --outdir        results/
 */
workflow {

    feature_table_ch = Channel.fromPath(params.feature_table)
    taxonomy_ch      = Channel.fromPath(params.taxonomy)
    metadata_ch      = Channel.fromPath(params.metadata)

    alpha_ch = params.containsKey('alpha') && params.alpha
        ? Channel.fromPath(params.alpha)
        : Channel.of(file('NO_FILE'))

    GROUPD(
        feature_table_ch,
        taxonomy_ch,
        metadata_ch,
        alpha_ch,
    )
    // Note: report_generator/ is passed inside GROUPD via projectDir
}
