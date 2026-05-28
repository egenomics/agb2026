process MULTIQC {
    container "containers/multiqc_1.19.sif"

    // Publishes the HTML and Data folder to the quality metrics directory
    publishDir "${params.out_dir}/quality/original_multiqc", mode: 'copy', pattern: "multiqc_report.html"
    publishDir "${params.out_dir}/quality/original_multiqc", mode: 'copy', pattern: "multiqc_data"

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data_dir

    // Isolate the text file here so we can pass it directly to the Python script
    path "multiqc_fastqc.txt",  emit: fastqc_txt

    script:
    """
    multiqc .

    # Extract the file to the top level of the working directory so Nextflow can grab it
    cp multiqc_data/multiqc_fastqc.txt multiqc_fastqc.txt
    """
}
