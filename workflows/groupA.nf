/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GROUP A: DATA HANDLING & PREPROCESSING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This workflow orchestrates all modules from Group A.
    Responsibilities:
    - Data input and validation
    - Quality control
    - Preprocessing and normalization
    - Metadata integration

    PLACEHOLDER: Add your group's modules and subworkflows below
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO: Import your group's modules here
// include { MODULE_NAME } from '../modules/groupA/module_name/main'
// include { SUBWORKFLOW_NAME } from '../subworkflows/local/groupA_subworkflow'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GROUPA {

    take:
    ch_samplesheet  // channel: samplesheet read in from --input
    ch_metadata     // channel: metadata file

    main:

    ch_versions = channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PLACEHOLDER: Add your team's workflow logic here
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Example structure (replace with your actual modules):
    //
    // QUALITY_CHECK(ch_samplesheet)
    // ch_versions = ch_versions.mix(QUALITY_CHECK.out.versions.first())
    //
    // PREPROCESSING(QUALITY_CHECK.out.filtered_data)
    // ch_versions = ch_versions.mix(PREPROCESSING.out.versions.first())
    //
    // METADATA_INTEGRATION(ch_metadata, QUALITY_CHECK.out.sample_ids)
    // ch_versions = ch_versions.mix(METADATA_INTEGRATION.out.versions.first())

    emit:

    // TODO: Define your team's outputs
    // trimmed_reads = PREPROCESSING.out.trimmed           // channel: Preprocessed sequences
    // metadata = METADATA_INTEGRATION.out.metadata         // channel: Integrated metadata
    // qc_reports = QUALITY_CHECK.out.reports              // channel: QC reports

    versions = ch_versions                                  // channel: Software versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
