process CLEAN_MULTIQC {
    // Load the needed containers
    // pandas: the galaxyproject depot has no modern build, but quay.io is reachable
    // from the cluster and singularity pulls it directly. pandas 1.5.2 is API-compatible
    // with the simple read_csv/to_csv used by clean_multiqc.py.
    container 'quay.io/biocontainers/pandas:1.5.2'
    // publishDir handled by the generic resolver in conf/modules.config -> ${outdir}/groupA/clean/
    input:
    // Receives the raw text summary file emitted by the MultiQC process
    path multiqc_fastqc_txt
    output:
    // Emits the final cleaned table
    path "quality_report_A.tsv", emit: final_report

    script:
    """
    ${projectDir}/bin/groupA/clean_multiqc.py ${multiqc_fastqc_txt} quality_report_A.tsv
    """
}
