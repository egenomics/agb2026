process CUTADAPT {
    // Defines the tags for nexftlow output in the terminal
    tag "${sample_id}"
    // Load the needed containers
    // Fallback (local image): containers/cutadapt_4.6.sif
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:4.6--py39hf95cd2a_1' :
        'biocontainers/cutadapt:4.6--py39hf95cd2a_1' }"
    input:
    tuple val(sample_id), path(reads)
    output:
    // Emits the trimmed reads (gzipped, as Group B expects) alongside the sample_id
    tuple val(sample_id), path("${sample_id}_trim.fastq.gz"), emit: trimmed_reads

    script:
    """
    cutadapt -g ${params.fwd_primer} -o ${sample_id}_trim.fastq.gz --cores ${task.cpus} ${reads[0]}
    """
}
