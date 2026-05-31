process METADATA_REPORT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path (metadata)

    output:
    path "clinical_association_map", emit: co_heat

    script:
    """
    python co_heat.py ${metadata} "clinical_association_map"
    """
}