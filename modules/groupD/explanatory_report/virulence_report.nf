
process VIRULENCE_PLOT {
    publishDir "${params.outdir}/explanatory", mode: 'copy'

    container "containers/groupD.sif"

    input:
    path meta 
    path asv_table 
    path annotated_taxonomy 
    

    output:
    path "Virulence_analysis", emit: virulence_report

    script:
    """
    virulence_report.sh ${asv_table} ${annotated_taxonomy} ${meta}
    """
}
