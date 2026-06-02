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
    // Publish the final compiled website to your results directory
    publishDir "${params.outdir}", mode: 'copy'

    label "groupD"

    input:
    // This catches an aggregated list of all plot directories
    path plot_dirs, stageAs: 'plot_dir*/*'

    output:
    // The outputs created by build_album.py
    path "report", emit: report

    script:
    """
    build_report.sh report $plot_dirs
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
