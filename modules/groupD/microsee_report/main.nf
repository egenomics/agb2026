process MICROSEE_REPORT {
    tag "microsee_report"
    label 'process_low'

    conda "conda-forge::python>=3.10,<3.13 conda-forge::pandas>=2.2,<3 conda-forge::numpy>=1.26,<3 conda-forge::pydantic>=2.7,<3"
    container 'ghcr.io/agb2026/microsee:latest'

    publishDir params.outdir, mode: 'copy'

    input:
    path report_generator, stageAs: 'report_generator'   // stages the whole package into the work dir
    path feature_table
    path taxonomy
    path metadata
    path alpha         // pass file('NO_FILE') if not available

    output:
    path "microsee_report.html", emit: report

    script:
    def alpha_arg = alpha.name != 'NO_FILE' ? "--alpha ${alpha}" : ""
    """
    python3 report_generator/generate_report.py \\
        --feature-table ${feature_table} \\
        --taxonomy      ${taxonomy}      \\
        --metadata      ${metadata}      \\
        ${alpha_arg}                     \\
        --output        microsee_report.html
    """
}
