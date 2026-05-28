process CLEAN_MULTIQC {
    // Uses a lightweight container that has Python 3 and Pandas pre-installed
    container "containers/pandas_2.0.3.sif"

    // Saves the final cleaned TSV file into your quality output directory
    publishDir "${params.out_dir}/quality", mode: 'copy'

    input:
    // Receives the raw text summary file emitted by the MultiQC process
    path multiqc_fastqc_txt

    output:
    // Emits the final cleaned table
    path "quality_report_A.tsv", emit: final_report

    script:
    """
    # Nextflow automatically finds 'clean_multiqc.py' because it is in your bin/ folder.
    # We pass the input file and the desired output name as arguments.
    a_clean_multiqc.py ${multiqc_fastqc_txt} quality_report_A.tsv
    """
}
