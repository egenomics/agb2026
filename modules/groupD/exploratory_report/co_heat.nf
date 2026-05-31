process METADATA_REPORT {
    publishDir "${params.outdir}/exploratory", mode: 'copy'

    container "containers/groupD.sif"

    input:
    tuple path (metadata)

    output:
    path "clinical_association_map", emit: co_heat

    script:
    """
    python co_heat.py ${metadata} "clinical_association_map"
    """
}