process CUTADAPT {
    tag "${sample_id}"

    // Points directly to your local Singularity container
    container "${projectDir}/containers/cutadapt_4.6.sif"

    // Saves trimmed fastq files directly to your sequence directory
    publishDir "${params.out_dir}/seqs/trimmed", mode: 'copy'

    input:
    // Takes the sample identifier and a list containing [read_1]
    tuple val(sample_id), path(reads)

    output:
    // Emits the trimmed pair alongside the sample_id for upstream processes like Bowtie2 or FastQC
    tuple val(sample_id), path("${sample_id}_trim.fastq"), emit: trimmed_reads

    script:
    """
    cutadapt \\
        -g ${params.fwd_primer} \\
        -o ${sample_id}_trim.fastq \\
        --cores ${task.cpus} \\
        ${reads[0]}
    """
}
