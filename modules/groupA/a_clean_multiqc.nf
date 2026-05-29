process CLEAN_MULTIQC {
    // Load the needed containers
    container "containers/pandas_2.0.3.sif"
    // Publish the results into the quality report directory
    publishDir "${params.out_dir}/quality_report", mode: 'copy'
    input:
    // Receives the raw text summary file emitted by the MultiQC process
    path multiqc_fastqc_txt
    output:
    // Emits the final cleaned table
    path "quality_report_A.tsv", emit: final_report

    script:
    """
    groupA/clean_multiqc.py ${multiqc_fastqc_txt} quality_report_A.tsv
    """
}
