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
    
    emit:
    alpha_plots = EXPLANATORY_REPORT.out.alpha_plots
    pcoa_plots  = EXPLANATORY_REPORT.out.pcoa_plots 

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
