/*
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

include { EXPLANATORY_REPORT }         from '../subworkflows/groupD_explanatory_report'
include { EXPLORATORY_REPORT }         from '../subworkflows/groupD_exploratory_report'
include { MICROSEE_REPORT }    from '../modules/groupD/reporting_module/report'
include { MICROSEE_PATIENT_REPORT }    from '../modules/groupD/reporting_module/report'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GROUPD {

    take:
    ch_metadata  // channel: Metadata from Group A
    ch_alpha_tuple    // channel: alpha_div Metrics
    ch_beta_tuple  // channel: beta_div Metrics
    ch_rarefaction // channel: the plots from rarefaction
    ch_annotated_counts // channel: Path to "annotated_table_counts.tsv"
    ch_annotated_taxonomy // channel: Path to "annotated_taxonomy.tsv"
    ch_summary // channel: Path to "contamination_summary.tsv"

    main:

    // Extract only the Bray Diversity
    ch_bray = ch_beta_tuple.map { it[0] }

    EXPLORATORY_REPORT(ch_metadata, ch_annotated_counts, ch_annotated_taxonomy, ch_bray)

    EXPLANATORY_REPORT(
        ch_metadata,
        ch_annotated_counts,                      // Counts for each sample
        ch_annotated_taxonomy,                    // Taxonomies in counts
        ch_alpha_tuple,                           // tuple with all alphas
        EXPLORATORY_REPORT.out.pca_object,        // PCA RDS path
        EXPLORATORY_REPORT.out.pca_genus_counts   // PCA genus counts
    )


    // ---------------------------------------------------------
    // 3. Prepare and report building process
    // ---------------------------------------------------------
    // This will need a collection of all the plot dirs and files
    ch_all_plots = channel.empty()
        .mix(
            EXPLORATORY_REPORT.out.pca_results,
            EXPLORATORY_REPORT.out.parallel_plot,
            EXPLORATORY_REPORT.out.heatmap_plot,
            EXPLORATORY_REPORT.out.volcano_plot,
            EXPLORATORY_REPORT.out.beta_tree,
            ch_rarefaction,
            ch_summary,
            EXPLANATORY_REPORT.out.overview_table,
            EXPLANATORY_REPORT.out.alpha_plots,
            EXPLANATORY_REPORT.out.pcoa_plots,
            EXPLANATORY_REPORT.out.patogeny_plots,
            EXPLANATORY_REPORT.out.violin_plots
        )
        .collect()

    MICROSEE_REPORT(ch_all_plots)
    
    emit:
    alpha_plots = EXPLANATORY_REPORT.out.alpha_plots
    pcoa_plots  = EXPLANATORY_REPORT.out.pcoa_plots
    html_report = MICROSEE_REPORT.out.report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
