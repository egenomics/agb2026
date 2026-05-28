process FASTQC {
    tag "${sample_id} (${stage})"

    // Points directly to your local Singularity container
    container "containers/fastqc_0.12.1.sif"

    input:
    tuple val(sample_id), val(stage), path(reads)

    output:
    // Emits the HTML and ZIP files for downstream MultiQC processing
    tuple val(sample_id), val(stage), path("*_fastqc.{html,zip}"), emit: qc_files

    script:
    """
    fastqc -o . --threads ${task.cpus} ${reads}
    """
}
