process VIOLIN_PLOTS {
    publishDir "${params.outdir}/explanatory", mode: 'copy'

    input:
    path meta 
    path asv_table 
    path annotated_taxonomy 
    

    output:
    path("violin_analysis"), emit: violin_report

    script:
    """
    violin.sh ${asv_table} ${annotated_taxonomy} ${meta} "violin_analysis"
    """
}