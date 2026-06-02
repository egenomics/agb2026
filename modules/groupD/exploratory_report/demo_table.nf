process DEMOGRAPHIC_SUMMARY {
    publishDir "${params.outdir}/exploratory", mode: 'copy'

    label "groupD"

    input:
    path(metadata)

    output:
    path "demographic_table", emit: demo_table_png

    script:
    """
    demo_table.py ${metadata} "demographic_table"
    """
}