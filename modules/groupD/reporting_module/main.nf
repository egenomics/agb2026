// ─────────────────────────────────────────────────────────────────────────────
// Group D — MicroSee report process
// AGB 2026
//
// Final "make it pretty" step: collects the plot directories produced by the
// exploratory and explanatory processes and assembles a single self-contained
// HTML report (index.html + one standalone page per patient).
//
// The generator (bin/build_album.py) scans recursively and routes:
//   - PNGs whose filename contains an ERR id  -> per-patient report
//   - all other PNGs                          -> cohort (exploratory) report
//   - non-PNG files (.rds, .tsv, .qza)        -> ignored
// So it does not matter how the upstream outputs are nested or named.
//
// bin/ is auto-staged onto PATH by Nextflow, matching the other Group D procs.
// ─────────────────────────────────────────────────────────────────────────────

process MICROSEE_REPORT {

    tag "microsee_report"
    publishDir "${params.outdir}", mode: 'copy'
    conda 'conda-forge::python=3.11'

    input:
    // Any number of plot directories from the upstream processes. Collect them
    // into one channel in the workflow and pass as a list, e.g.:
    //   MICROSEE_REPORT( plots_ch.collect() )
    path plot_dirs

    output:
    path "report", emit: report

    script:
    """
    build_report.sh report ${plot_dirs}
    echo "✓ MicroSee report at \${PWD}/report/index.html"
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// Optional per-patient process — builds ONE patient's standalone report.
// ─────────────────────────────────────────────────────────────────────────────
process MICROSEE_PATIENT_REPORT {

    tag { patient_id }
    publishDir "${params.outdir}/report/patients", mode: 'copy'
    conda 'conda-forge::python=3.11'

    input:
    tuple val(patient_id), path(plot_dirs)

    output:
    path "${patient_id}.html", emit: patient_report

    script:
    """
    build_report.sh . ${plot_dirs} -- ${patient_id}
    echo "✓ Patient report: ${patient_id}.html"
    """
}
