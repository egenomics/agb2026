process NEW_SAMPLE_SHEET {
    // Load container
    // singularity pulls it directly. 1.5.2 is fine for generate_sample_sheet_B.py.
    container 'quay.io/biocontainers/pandas:1.5.2'
    input:
    path id_list

    output:
    path "sample_sheet_B.csv", emit: final_sheet

    script:
    """
    # Call the script from the bin/ folder.
    # Argument 1: The input text file (${id_list})
    # Argument 2: The dynamic path where the trimmed files are saved
    # Argument 3: The name of the output CSV
    
    ${projectDir}/bin/groupA/generate_sample_sheet_B.py ${id_list} "${params.outdir}/seqs/trimmed" sample_sheet_B.csv
    """
}