process PARALLEL_PLOT_MICROBIOME {

    tag "${metadata.baseName}"

    publishDir "${params.outdir}/exploratory",
        mode: 'copy'

    input:
    path(asv_table)
    path(taxonomy)
    path(metadata)

    output:
    path("parallel_results"), emit: parallel_results

    script:
    """
    mkdir -p parallel_results

    Rscript parallel_plot_relabu_healthygrouped.R \
        --asv ${asv_table} \
        --taxonomy ${taxonomy} \
        --metadata ${metadata} \
        --outdir parallel_results
    """
}