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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GROUPD {

    take:
    ch_metadata  // channel: Metadata from Group A
    ch_groupC    // channel: Validated results from Group C
    ch_groupB  // channel: Analysis data

    main:

    ch_versions = channel.empty()

    //Assuming groupA sends metadata and groupC the diversities
    ch_input_for_alpha = ch_metadata.combine(ch_groupC)

    // The channel is just passed to the process (it unwraps the tuple)
    alpha_diversity_status(ch_input_for_alpha)

    emit:
    // Captures the process `.out` emitted as the name `alpha_div_dist_dir`
    alpha_plots = alpha_diversity_status.out.alpha_div_dist_dir

    versions = ch_versions // channel: Software versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
