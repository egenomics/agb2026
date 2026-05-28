process CUTADAPT {
    // Defines the tags for nexftlow output in the terminal
    tag "${sample_id}"
    // Load the needed containers
    container "/containers/cutadapt_4.6.sif"
    input:
    tuple val(sample_id), path(reads)
    output:
    // Emits the trimmed pair alongside the sample_id
    tuple val(sample_id), path("${sample_id}_trim.fastq"), emit: trimmed_reads

    script:
    """
    cutadapt -g ${params.fwd_primer} -o ${sample_id}_trim.fastq --cores ${task.cpus} ${reads[0]}
    """
}
