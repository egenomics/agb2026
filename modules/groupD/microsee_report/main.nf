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
    path alpha   // pass file('NO_FILE') if not available
    val  mode    // 'cohort' (default) | 'patient' | 'all'

    output:
    path "microsee_report.html",  emit: report                        // primary cohort report
    path "microsee_report_*.html", emit: patient_reports, optional: true  // per-patient (mode=patient|all)

    script:
    // We invoke generate_report.py directly rather than the microsee-report CLI entrypoint
    // because the CLI is only available after `pip install -e .`, which is not run inside the
    // Nextflow work directory. python3 on the staged script is equivalent and works in all envs.
    def alpha_arg = alpha.name != 'NO_FILE' ? "--alpha ${alpha}" : ""
    def mode_arg  = mode ? "--mode ${mode}" : "--mode cohort"
    """
    # Verify that the bundled Plotly.js is present — without it all charts will be blank.
    # This file must be committed to git; it cannot be downloaded on offline HPC nodes.
    if [ ! -f report_generator/charts/plotly.min.js ]; then
        echo "ERROR: report_generator/charts/plotly.min.js not found." >&2
        echo "Commit it to git first — see README Troubleshooting section." >&2
        echo "  curl -fsSL https://cdn.plot.ly/plotly-2.35.2.min.js \\" >&2
        echo "       -o modules/groupD/microsee_report/report_generator/charts/plotly.min.js" >&2
        echo "  git add modules/groupD/microsee_report/report_generator/charts/plotly.min.js" >&2
        echo "  git commit -m 'Bundle Plotly.js v2.35.2 for offline HPC use'" >&2
        exit 1
    fi

    python3 report_generator/generate_report.py \\
        --feature-table ${feature_table} \\
        --taxonomy      ${taxonomy}      \\
        --metadata      ${metadata}      \\
        ${alpha_arg}                     \\
        ${mode_arg}                      \\
        --output        microsee_report.html
    """
}
