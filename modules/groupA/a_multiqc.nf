process MULTIQC {
    // Load the needed containers
    container "containers/multiqc_1.19.sif"

    // Publishes the results into the quality report directory
    publishDir "${params.out_dir}/quality_report", mode: 'copy', pattern: "multiqc_report.html"
    publishDir "${params.out_dir}/quality_report", mode: 'copy', pattern: "multiqc_data"

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data_dir
    path "multiqc_fastqc.txt",  emit: fastqc_txt

    script:
    """
    multiqc .
    """
}
