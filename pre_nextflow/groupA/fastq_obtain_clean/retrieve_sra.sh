#!/bin/bash

# Prepare a conda enviorment to run the function (or your own alternative)
# conda create -n retrieve_sra_env -y
conda activate retrieve_sra_env
# conda install bioconda::entrez-direct -y

# 1. Help function
show_help() {
    cat << EOF
Usage: ${0##*/} [-h] -i INPUT_FILE [-o OUTPUT_FILE]

Description: Reads a metadata TSV file, searches NCBI SRA using the 'sample_name_id' column as a full string query, and prepends a new 'sample-id' column containing the resulting SRA IDs.
It also identifies, reports, and removes any samples that could not be resolved (NA).

OPTIONS:
    -h         Show this help message.
    -i FILE    (Required) Path to the input metadata file (TSV).
    -o FILE    (Optional) Path to the output metadata file. Defaults to appending '_clean' to the input filename.
EOF
}

# 2. Parse the argument into the script
INPUT_FILE=""
OUTPUT_FILE=""

# The colon after 'i' and 'o' means they require an argument
while getopts "hi:o:" opt; do
    case "$opt" in
        h) show_help; exit 0 ;;
        i) INPUT_FILE=$OPTARG ;;
        o) OUTPUT_FILE=$OPTARG ;;
        *) show_help; exit 1 ;;
    esac
done

# 3. Ensure the required input file is provided
if [ -z "$INPUT_FILE" ]; then
    echo "Error: You must provide an input file using the -i flag."
    echo ""
    show_help
    exit 1
fi

# 4. Check if input file actually exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found."
    exit 1
fi

# 4b. Set the output name if not provided by the user
if [ -z "$OUTPUT_FILE" ]; then
    OUTPUT_FILE="${INPUT_FILE%.*}_clean.${INPUT_FILE##*.}"
fi

echo "Processing '$INPUT_FILE'..."

# 5. Find the 'sample_name_id' column
# Hardcode the delimiter for TSV (tab)
DELIM=$'\t'
TMP_SRA="sra_temp_column.txt"

# Extract header and strip hidden Windows carriage returns (\r) and quotes
header=$(head -n 1 "$INPUT_FILE" | tr -d '\r' | tr -d '"')
col_idx=$(echo "$header" | tr "$DELIM" '\n' | grep -nx "sample_name_id" | cut -d: -f1)

# Safety check
if [ -z "$col_idx" ]; then
    echo "Error: Could not find 'sample_name_id' column in the header."
    exit 1
fi

# Write the header for our temporary column
echo "sample-id" > "$TMP_SRA"

# 6. Loop and fetch the sra code ids
tail -n +2 "$INPUT_FILE" | cut -d"$DELIM" -f"$col_idx" | while read -r sample_raw; do
    
    # Clean the sample string: Remove \r, quotes, and extra spaces
    sample=$(echo "$sample_raw" | tr -d '\r' | tr -d '"' | xargs)

    # Skip empty lines
    if [ -z "$sample" ]; then
        echo "NA" >> "$TMP_SRA"
        continue
    fi
    
    echo "  -> Fetching: $sample"
    
    # Make sure this works properly
    sra_id=$(esearch -db sra -query "$sample" </dev/null | efetch -format runinfo 2>/dev/null | cut -d ',' -f 1 | grep -E '^[SED]RR[0-9]+' | head -n 1)

    # If the string is empty, put NA, otherwise put the SRA ID
    if [ -z "$sra_id" ]; then
        echo "NA" >> "$TMP_SRA"
    else
        echo "$sra_id" >> "$TMP_SRA"
    fi
    
    # Small pause to avoid overloading NCBI servers
    sleep 1

done

# 7. Merge data, placing the new column as the first one
echo "Merging data..."
# By putting TMP_SRA first in the paste command, 'sample-id' becomes the first column
paste -d"$DELIM" "$TMP_SRA" "$INPUT_FILE" > "$OUTPUT_FILE"

# Clean up
rm "$TMP_SRA"

# 8. Handle NA rows: Display, log, and eliminate them
echo "Checking for unresolved (NA) samples..."

FAILED_LOG="${OUTPUT_FILE%.*}_failed_samples.txt"
# Since we added 1 column at the beginning, the sample_name_id column shifted right by 1
NEW_COL_IDX=$((col_idx + 1))

# Extract 'sample_name_id' for rows where the first column is "NA"
awk -F"$DELIM" -v c="$NEW_COL_IDX" '$1 == "NA" {print $c}' "$OUTPUT_FILE" > "$FAILED_LOG"

# Check if the failed log has any contents
if [ -s "$FAILED_LOG" ]; then
    echo "---------------------------------------------------"
    echo "Warning: The following sample_name_ids returned 'NA':"
    cat "$FAILED_LOG" | while read -r failed_name; do
        echo "  - $failed_name"
    done
    echo "---------------------------------------------------"
    echo "These names have been saved to: $FAILED_LOG"
    
    # Remove the NA rows from the output file using awk, output to a temp file, then overwrite
    echo "Removing NA rows from the final output..."
    awk -F"$DELIM" '$1 != "NA"' "$OUTPUT_FILE" > "${OUTPUT_FILE}.tmp"
    mv "${OUTPUT_FILE}.tmp" "$OUTPUT_FILE"
else
    echo "All samples successfully resolved."
    # Remove the empty log file since it's not needed
    rm -f "$FAILED_LOG"
fi

echo "Finished! Cleaned file saved to: $OUTPUT_FILE"
