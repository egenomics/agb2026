process kraken2 {
    tag "kraken2_classification"

    input:
    path rep_seqs

    output:
    path "kraken2_report.txt",     emit: report
    path "kraken2_raw_output.txt", emit: raw_output

    script:
    """
    kraken2 --db ${params.kraken2_db} \\
            --threads ${task.cpus} \\
            --use-names \\
            --report kraken2_report.txt \\
            ${rep_seqs} > kraken2_raw_output.txt
    """
}

workflow {
    rep_seqs = file("${params.data_dir}/rep_seqs.fasta")

    kraken2(rep_seqs)
}