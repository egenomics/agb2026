process VIRULENCE_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path meta 
    path asv_table 
    path annotated_taxonomy 
    

    output:
    path "Violin Analysis", emit: violin_report

    script:
    """
    violin.sh ${asv_table} ${annotated_taxonomy} ${meta}
    ""