process MICROSEE_REPORT {
    tag "microsee_report"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.3.3 conda-forge::numpy=2.4.1 conda-forge::pydantic=2.13.4"
    container 'ghcr.io/egenomics/microsee-report:latest'

    publishDir params.outdir, mode: 'copy'

    input:
    path report_generator, stageAs: 'report_generator'   // stages the whole package into the work dir
    path feature_table
    path taxonomy
    path metadata
    path alpha            // pass file('NO_FILE') if not available
    path distance_matrix  // pass file('NO_FILE') if not available
    val  mode             // 'cohort' (default) | 'patient' | 'all'

    output:
    path "microsee_report.html",  emit: report                        // primary cohort report
    path "microsee_report_*.html", emit: patient_reports, optional: true  // per-patient (mode=patient|all)

    script:
    // We invoke generate_report.py directly rather than the microsee-report CLI entrypoint
    // because the CLI is only available after `pip install -e .`, which is not run inside the
    // Nextflow work directory. python3 on the staged script is equivalent and works in all envs.
    def alpha_arg = alpha.name != 'NO_FILE' ? "--alpha ${alpha}" : ""
    def dist_arg  = distance_matrix.name != 'NO_FILE' ? "--distance-matrix ${distance_matrix}" : ""
    def mode_arg  = mode ? "--mode ${mode}" : "--mode cohort"
    """
    # Verify that the bundled Plotly.js is present — without it all charts will be blank.
    # This file must be committed to git; it cannot be downloaded on offline HPC nodes.
    if [ ! -f report_generator/charts/plotly.min.js ]; then
        echo "ERROR: report_generator/charts/plotly.min.js not found." >&2
        echo "This file must be committed to git — see README Troubleshooting." >&2
        exit 1
    fi

    python3 report_generator/generate_report.py \\
        --feature-table ${feature_table} \\
        --taxonomy      ${taxonomy}      \\
        --metadata      ${metadata}      \\
        ${alpha_arg}                     \\
        ${dist_arg}                      \\
        ${mode_arg}                      \\
        --output        microsee_report.html
    """
}
