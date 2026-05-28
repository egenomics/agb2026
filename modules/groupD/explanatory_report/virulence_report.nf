
process VIRULENCE_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path metadata 
    path asv_table 
    path annotated_taxonomy 
    

    output:
    path "Virulence_analysis", emit: virulence_report

    script:
    """
    virulence_report.sh ${asv_table} ${annotated_taxonomy} ${metadata}
    """
}
