process contamination_filter {
    tag "contamination_filtering"

    input:
    path asv_table
    path taxonomy
    path kraken2_output
    path blanks_metadata
    path patients_metadata

    output:
    path "annotated_table_counts.tsv", emit: annotated_counts
    path "annotated_taxonomy.tsv",     emit: annotated_taxonomy
    path "contamination_summary.tsv",  emit: summary

    script:
    """
    Rscript ${projectDir}/contamination_filtering.R \\
        ${asv_table} \\
        ${taxonomy} \\
        ${kraken2_output} \\
        ${blanks_metadata} \\
        ${patients_metadata}
    """
}

workflow {
    asv_table         = file("${params.data_dir}/asv_table.tsv")
    taxonomy          = file("${params.data_dir}/ASV_taxonomy.tsv")
    kraken2_output    = file("${params.data_dir}/kraken2_raw_output.txt")
    blanks_metadata   = file("${params.data_dir}/blank_sample_information.tsv")
    patients_metadata = file("${params.data_dir}/patients_sample_information.tsv")

    contamination_filter(asv_table, taxonomy, kraken2_output, blanks_metadata, patients_metadata)
}