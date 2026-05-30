process MULTIQC {
    // Load the needed containers
    // Fallback (local image): containers/multiqc_1.19.sif
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.19--pyhdfd78af_0' :
        'biocontainers/multiqc:1.19--pyhdfd78af_0' }"

    // publishDir handled by the generic resolver in conf/modules.config -> ${outdir}/groupA/multiqc/

    input:
    path qc_files

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data_dir
    path "multiqc_data/multiqc_fastqc.txt",  emit: fastqc_txt

    script:
    """
    multiqc .
    """
}
