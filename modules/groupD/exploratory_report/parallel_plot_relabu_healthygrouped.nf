process PARALLEL_PLOT_MICROBIOME {

    tag "${metadata.baseName}"

    publishDir "${params.outdir}/parallel_plot",
        mode: 'copy'

    container "containers/r_groupD.sif"

    input:
    path(asv_table)
    path(taxonomy)
    path(metadata)

    output:
    path("parallel_results")

    script:
    """
    mkdir -p parallel_results

    Rscript ${moduleDir}/parallel_plot_relabu_healthygrouped.R \
        --asv ${asv_table} \
        --taxonomy ${taxonomy} \
        --metadata ${metadata} \
        --outdir parallel_results
    """
}