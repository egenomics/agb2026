process VOLCANO_PLOT {
    publishDir "${params.outdir}/exploratory", mode: 'copy'

    input:
    path(metadata)  
    path(asv_table)
    path(asv_taxonomy)
    

    output:
    path "volcano_plot", emit: volcano_png

    script:
    """
    volcano.R ${asv_table} ${asv_taxonomy} ${metadata} "volcano_plot"
    """
}