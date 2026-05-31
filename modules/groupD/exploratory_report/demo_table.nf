process DEMOGRAPHIC_SUMMARY {
    publishDir "${params.outdir}/exploratory", mode: 'copy'

    container "containers/groupD.sif"

    input:
    tuple path metadata

    output:
    path "demographic_table", emit: demo_table_png

    script:
    """
    python demo_table.py ${metadata} "demographic_table"
    """
}