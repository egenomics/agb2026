process METADATA_REPORT {
    publishDir "${params.outdir}/exploratory", mode: 'copy'

    label "groupD"

    input:
    path(metadata)

    output:
    path "clinical_association_map", emit: co_heat

    script:
    """
    co_heat.py ${metadata} "clinical_association_map"
    """
}