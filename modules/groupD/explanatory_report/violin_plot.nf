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
    violin.sh --asv-table ${asv_table} --taxonomy ${annotated_taxonomy} --metadata ${meta} --outdir .
    """
}