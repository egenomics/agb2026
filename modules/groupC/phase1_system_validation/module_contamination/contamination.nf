process KRAKEN2 {
    tag "kraken2_classification"
    label 'process_medium'

    conda "bioconda::kraken2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kraken2:2.1.3--pl5321hdcf5f25_0' :
        'biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0' }"

    input:
    path rep_seqs

    output:
    path "kraken2_report.txt",     emit: report
    path "kraken2_raw_output.txt", emit: raw_output

    script:
    def args = task.ext.args ?: ''
    """
    kraken2 --db ${params.kraken2_db} \\
            --threads ${task.cpus} \\
            --use-names \\
            --report kraken2_report.txt \\
            ${args} \\
            ${rep_seqs} > kraken2_raw_output.txt
    """
}

process CONTAMINATION_FILTER {
    tag "contamination_filtering"
    label 'process_low'

    conda "bioconda::bioconductor-decontam conda-forge::r-tidyverse conda-forge::r-gridextra"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-decontam_r-gridextra_r-tidyverse:3cd9a275a53c9075' :
        'community.wave.seqera.io/library/bioconductor-decontam_r-gridextra_r-tidyverse:3cd9a275a53c9075' }"

    input:
    path asv_table
    path taxonomy
    path kraken2_raw_output
    path metadata

    output:
    path "contamination_summary.png", emit: summary_png

    script:
    def args = task.ext.args ?: ''
    """
    contamination_filtering.R \\
        ${asv_table} \\
        ${taxonomy} \\
        ${kraken2_raw_output} \\
        ${metadata}
    """
}

workflow CONTAMINATION {

    take:
    asv_table  // path: asv_table.tsv from Group B
    taxonomy   // path: ASV_taxonomy.tsv from Group B
    rep_seqs   // path: rep_seqs.fasta from Group B
    metadata   // path: sample-metadata.tsv from Group A

    main:
    KRAKEN2(rep_seqs)

    CONTAMINATION_FILTER(
        asv_table,
        taxonomy,
        KRAKEN2.out.raw_output,
        metadata
    )

    emit:
    summary_png = CONTAMINATION_FILTER.out.summary_png  // contamination_summary.png
}