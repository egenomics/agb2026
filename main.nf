#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/abgtemplate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/abgtemplate
    Website: https://nf-co.re/abgtemplate
    Slack  : https://nfcore.slack.com/channels/abgtemplate
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// GROUP WORKFLOWS
include { GROUPA }  from './workflows/groupA'
include { GROUPB }  from './workflows/groupB'
include { GROUPC }  from './workflows/groupC'
include { GROUPD }  from './workflows/groupD'

// PIPELINE UTILITIES
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_abgtemplate_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_abgtemplate_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Orchestrate all group workflows
//
workflow NFCORE_ABGTEMPLATE {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    ch_metadata = file(params.metadata, checkIfExists: true)

    //
    // GROUP A: Data Handling & Preprocessing
    //
    GROUPA(
        samplesheet,
        ch_metadata
    )

    //
    // GROUP B: Analysis
    //
    ch_multiqc_config = file("${projectDir}/assets/multiqc_config.yml", checkIfExists: false)
    ch_multiqc_logo   = file("${projectDir}/assets/multiqc_logo.png", checkIfExists: false)

    GROUPB(
        GROUPA.out.trimmed_reads,          // [meta, [reads]] — per-sample handoff from Group A
        ch_multiqc_config,
        ch_multiqc_logo,
        channel.empty(),                   // ch_collated_versions   (QC_CHECKS disabled in Group B)
        channel.empty(),                   // ch_methods_description
        channel.empty(),                   // ch_workflow_summary
        file(params.dada2_train_set),
        file(params.dada2_species_set)
    )

    //
    // GROUP C / D: Validation & Reporting — NOT yet wired.
    // Their workflows exist but currently emit only `versions`; the channels
    // main.nf would consume (analysis_results, validated_results, html_report)
    // are not produced yet. Re-enable once Groups C and D implement their outputs.
    //
    // GROUPC(GROUPB.out.table_counts, GROUPA.out.metadata)
    // GROUPD(GROUPC.out.validated_results, GROUPB.out.table_counts)

    emit:
    multiqc_report = GROUPA.out.quality_report // channel: Group A MultiQC (until Group D produces the final report)
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_ABGTEMPLATE (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_ABGTEMPLATE.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
