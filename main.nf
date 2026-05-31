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
    GROUPB(
        GROUPA.out.trimmed_reads,
        GROUPA.out.metadata
    )

    //
    // GROUP C: Validation
    //
    GROUPC(
        GROUPB.out.analysis_results,
        GROUPA.out.metadata
    )

    //
    // GROUP D: Reporting
    //
        GROUPD(
            GROUPA.out.metadata, // channel: Metadata from Group A
            GROUPC.out.OutputMetricResultsAlpha,
            GROUPC.out.OutputMetricResultsBeta,
            GROUPC.out.OutputMetricResultsRarefaction,
            GROUPC.out.annotated_counts,
            GROUPC.out.annotated_taxonomy,
            GROUPC.out.summary /// channel: Path to "contamination_summary.tsv"
        )

    emit:
    multiqc_report = GROUPD.out.html_report // channel: /path/to/report
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
