process HEATMAP_MICROBIOME {

    tag "${metadata.baseName}"

    publishDir "${params.outdir}/exploratory",
        mode: 'copy'

    input:
    path(asv_table)
    path(taxonomy)
    path(metadata)

    output:
    path("heatmap_results"), emit: heatmap_results

    script:
    """
    mkdir -p heatmap_results

    heatmap_samples_bacteria.R \
        --asv ${asv_table} \
        --taxonomy ${taxonomy} \
        --metadata ${metadata} \
        --outdir heatmap_results
    """
}