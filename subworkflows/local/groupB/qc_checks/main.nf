// subworkflows/local/groupB/qc_checks/main.nf
//
// QC subworkflow: FastQC per-sample, aggregated by MultiQC.

include { FASTQC  } from '../../../../modules/nf-core/fastqc/main'
include { MULTIQC } from '../../../../modules/nf-core/multiqc/main'

workflow QC_CHECKS {

    take:
    ch_samplesheet             // channel: [ val(meta), [ reads ] ]
    multiqc_config             // path:    multiqc_config.yml
    multiqc_logo               // path:    multiqc_logo.png
    ch_collated_versions       // channel: [ path(versions.yml) ]
    ch_methods_description     // channel: [ path(methods_description_mqc.yaml) ]
    ch_workflow_summary        // channel: [ path(workflow_summary_mqc.yaml) ]

    main:
    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { _meta, file -> file })

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_collated_versions,
        ch_workflow_summary,
        ch_methods_description
    )

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [[id: 'alphaflow_qc'], files, multiqc_config, multiqc_logo, [], []]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report
    fastqc_zip     = FASTQC.out.zip
    versions       = ch_versions
}
