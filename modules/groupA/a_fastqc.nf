process FASTQC {
    // Defines the tags for nexftlow output in the terminal
    tag "${sample_id} (${stage})"
    // Load the needed containers
    // Fallback (local image): containers/fastqc_0.12.1.sif
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(sample_id), val(stage), path(reads)

    output:
    // Emits the HTML and ZIP files that will be used for downstream MultiQC processing
    tuple val(sample_id), val(stage), path("*_fastqc.{html,zip}"), emit: qc_files

    script:
    """
    fastqc -o . --threads ${task.cpus} ${reads}
    """
}
