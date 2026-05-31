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
include { alpha_diversity_status } from '../modules/groupD/explanatory_report/alpha_diversity_status'
include { PCoA_plots }             from '../modules/groupD/explanatory_report/PCoA' 
include { overview_table }         from '../modules/groupD/explanatory_report/overview_table.nf' 


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

    ch_versions = channel.empty()

    // ---------------------------------------------------------
    // 1. Prepare and Run Alpha Diversity
    // ---------------------------------------------------------
    // Combined with metadata
    ch_input_for_alpha = ch_metadata.combine(ch_alpha_tuple)
    
    alpha_diversity_status(ch_input_for_alpha)

    // ---------------------------------------------------------
    // 2. Prepare and Run PCoA
    // ---------------------------------------------------------
    // PCoA process expects: tuple path(pca), path(metadata)
    //Change for the PCA
    ch_input_for_pcoa = ch_annotated_counts.combine(ch_metadata)
    
    PCoA_plots(ch_input_for_pcoa)

    // ---------------------------------------------------------
    // 3. Prepare and run the overview table
    // ---------------------------------------------------------
    // overview process expects: tuple path(metadata), path(pca), path(z_scores), path(genus_counts)
    //Change for the PCA and annotated genus
    ch_overview = ch_metadata
        .join(ch_annotated_counts)
        .join(ch_annotated_taxonomy)
        .combine(ch_input_for_alpha.out.PCoA_patient_plots_dir)

    overview_table(ch_overview)

    emit:
    alpha_plots = alpha_diversity_status.out.alpha_div_dist_dir
    pcoa_plots  = PCoA_plots.out.PCoA_patient_plots_dir 
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
