
process VIRULENCE_PLOT {
    publishDir "${params.outdir}/explanatory", mode: 'copy'

    label "groupD"

    input:
    path meta 
    path asv_table 
    path annotated_taxonomy 
    

    output:
    path("Virulence_analysis"), emit: virulence_report

    script:
    """
    virulence_report.sh ${asv_table} ${annotated_taxonomy} ${meta}
    """
}
