process NEW_SAMPLE_SHEET {
    // Load the needed containers
    container "containers/pandas_2.0.3.sif"
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
    
    groupA/generate_sample_sheet_B.py ${id_list} "${params.out_dir}/seqs/trimmed" sample_sheet_B.csv
    """
}