process FASTQC {
    // Defines the tags for nexftlow output in the terminal
    tag "${sample_id} (${stage})"
    // Load the needed containers
    container "containers/fastqc_0.12.1.sif"

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
