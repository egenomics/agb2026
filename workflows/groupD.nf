/*
<<<<<<< HEAD
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
 *   params.alpha              path to alpha-diversity.tsv  (optional)
 *   params.distance_matrix    path to distance-matrix.tsv  (optional, e.g. Bray-Curtis)
 *   params.mode               'cohort' (default) | 'patient' | 'all'
 *   params.outdir         where to publish the HTML report
 */

include { MICROSEE_REPORT } from '../modules/groupD/microsee_report/main'
=======
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GROUP D: REPORTING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This workflow orchestrates all modules from Group D.
    Responsibilities:
    - Visualization
    - Report generation
    - Dashboard creation
    - Output formatting

    PLACEHOLDER: Add your group's modules and subworkflows below
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO: Import your group's modules here
// include { MODULE_NAME } from '../modules/groupD/module_name/main'
// include { SUBWORKFLOW_NAME } from '../subworkflows/local/groupD_subworkflow'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
>>>>>>> Group_A/main

workflow GROUPD {

    take:
<<<<<<< HEAD
    feature_table   // channel: path
    taxonomy        // channel: path
    metadata        // channel: path
    alpha             // channel: path  (can be file('NO_FILE'))
    distance_matrix   // channel: path  (can be file('NO_FILE'))
    mode              // val: 'cohort' | 'patient' | 'all'

    main:
    MICROSEE_REPORT(
        Channel.of(file("${moduleDir}/../modules/groupD/microsee_report/report_generator")),
        feature_table,
        taxonomy,
        metadata,
        alpha,
        distance_matrix,
        mode,
    )

    emit:
    report          = MICROSEE_REPORT.out.report
    patient_reports = MICROSEE_REPORT.out.patient_reports
}

/*
 * Standalone entrypoint — run this workflow directly:
 *
 *   nextflow run workflows/groupD.nf \
 *       --feature_table path/to/feature-table.tsv \
 *       --taxonomy      path/to/taxonomy.tsv       \
 *       --metadata      path/to/metadata.tsv       \
 *       [--alpha            path/to/alpha-diversity.tsv] \
 *       [--distance_matrix  path/to/distance-matrix.tsv] \
 *       [--mode             all]                         \
 *       --outdir            results/
 */
workflow {

    // Validate required parameters before attempting to stage any files.
    // Channel.fromPath on a null param gives a cryptic Nextflow error; this gives a clear one.
    if (!params.feature_table) error "Missing required parameter: --feature_table"
    if (!params.taxonomy)      error "Missing required parameter: --taxonomy"
    if (!params.metadata)      error "Missing required parameter: --metadata"

    feature_table_ch = Channel.fromPath(params.feature_table, checkIfExists: true)
    taxonomy_ch      = Channel.fromPath(params.taxonomy,      checkIfExists: true)
    metadata_ch      = Channel.fromPath(params.metadata,      checkIfExists: true)

    alpha_ch = params.containsKey('alpha') && params.alpha
        ? Channel.fromPath(params.alpha, checkIfExists: true)
        : Channel.of(file('NO_FILE'))

    distance_matrix_ch = params.containsKey('distance_matrix') && params.distance_matrix
        ? Channel.fromPath(params.distance_matrix, checkIfExists: true)
        : Channel.of(file('NO_FILE'))

    GROUPD(
        feature_table_ch,
        taxonomy_ch,
        metadata_ch,
        alpha_ch,
        distance_matrix_ch,
        params.mode,
    )
}
=======
    ch_validated_results    // channel: Validated results from Group C
    ch_analysis_data        // channel: Analysis data (if needed for visualization)

    main:

    ch_versions = channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PLACEHOLDER: Add your team's workflow logic here
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Example structure (replace with your actual modules):
    //
    // VISUALIZATION_MODULE(ch_validated_results)
    // ch_versions = ch_versions.mix(VISUALIZATION_MODULE.out.versions.first())
    //
    // REPORT_GENERATION(VISUALIZATION_MODULE.out.plots, ch_validated_results)
    // ch_versions = ch_versions.mix(REPORT_GENERATION.out.versions.first())
    //
    // DASHBOARD_MODULE(ch_validated_results, VISUALIZATION_MODULE.out.plots)
    // ch_versions = ch_versions.mix(DASHBOARD_MODULE.out.versions.first())

    emit:

    // TODO: Define your team's final outputs
    // html_report = REPORT_GENERATION.out.html             // channel: HTML report
    // plots = VISUALIZATION_MODULE.out.plots               // channel: Visualization files
    // dashboard = DASHBOARD_MODULE.out.dashboard           // channel: Interactive dashboard

    versions = ch_versions                                   // channel: Software versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
>>>>>>> Group_A/main
