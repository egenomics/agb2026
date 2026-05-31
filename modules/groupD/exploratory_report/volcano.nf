process VOLCANO_PLOT {
    publishDir "${params.outdir}/exploratory", mode: 'copy'

    container "containers/groupD.sif"

    input:
    tuple path (metadata)  
    path asv_table
    path asv_taxonomy
    

    output:
    path "volcano_plot", emit: volcano_png

    script:
    """
    Rscript volcano.R ${asv_table} ${asv_taxonomy} ${metadata} "volcano_plot"
    """
}