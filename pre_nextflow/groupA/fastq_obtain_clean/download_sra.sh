#!/bin/bash

# Prepare a conda enviorment to run the function (or your own alternative)
# conda create -n download_sra_env -y
conda activate download_sra_env
# conda install bioconda::sra-tools -y

# 1. Helo function
show_help() {
    cat << EOF
Usage: ${0##*/} [-h] -s SAMPLE_SHEET -o OUTPUT_DIR

Description: Dynamically parses a metadata TSV file for the 'sample-id' column,
             downloads the sequence runs using prefetch via Conda safely, 
             and flattens the output into the designated directory.

OPTIONS:
    -h         Show this help message.
    -s FILE    Path to the input metadata file (TSV format).
    -o DIR     Path to the final output directory for the .sra files.
EOF
}

# 2. Parse Arguments
SAMPLE_SHEET=""
OUT_DIR=""

while getopts "hs:o:" opt; do
    case "$opt" in
        h) show_help; exit 0 ;;
        s) SAMPLE_SHEET=$OPTARG ;;
        o) OUT_DIR=$OPTARG ;;
        *) show_help; exit 1 ;;
    esac
done

# Ensure mandatory arguments are provided
if [ -z "$SAMPLE_SHEET" ] || [ -z "$OUT_DIR" ]; then
    echo "Error: Missing required arguments (-s and -o are mandatory)."
    show_help
    exit 1
fi

# Ensure the sample sheet exists
if [ ! -f "$SAMPLE_SHEET" ]; then
    echo "Error: Sample sheet file '$SAMPLE_SHEET' not found."
    exit 1
fi

# 3. Setup Directories
mkdir -p "$OUT_DIR"
TMP_DOWNLOAD_DIR="${OUT_DIR}/prefetch_tmp"
mkdir -p "$TMP_DOWNLOAD_DIR"


# 4. Find the 'sample-id' Column Dynamically
DELIM=$'\t'
header=$(head -n 1 "$SAMPLE_SHEET" | tr -d '\r' | tr -d '"')
col_idx=$(echo "$header" | tr "$DELIM" '\n' | grep -nx "sample-id" | cut -d: -f1)

if [ -z "$col_idx" ]; then
    echo "Error: Could not find 'sample-id' column in the metadata header."
    exit 1
fi

echo "[INFO] Successfully mapped 'sample-id' to column index: $col_idx"

# 5. Iterate Over the TSV and Download
# tail -n +2 skips the header row entirely
tail -n +2 "$SAMPLE_SHEET" | cut -d"$DELIM" -f"$col_idx" | while read -r sample_raw; do
    
    # Clean hidden Windows carriage returns (\r) or wrapping quotes
    sample_id=$(echo "$sample_raw" | tr -d '\r' | tr -d '"' | xargs)
    
    # Skip empty lines
    if [ -z "$sample_id" ]; then
        continue
    fi

    echo "=========================================================="
    echo "[INFO] Processing Download for Run: $sample_id"
    echo "=========================================================="


    if prefetch "$sample_id" --output-directory "$TMP_DOWNLOAD_DIR" </dev/null; then
        
        EXPECTED_SRA="${TMP_DOWNLOAD_DIR}/${sample_id}/${sample_id}.sra"

        # Verify download, move it up to flatten folder tree, and remove subfolders
        if [ -f "$EXPECTED_SRA" ]; then
            echo "[INFO] Moving ${sample_id}.sra to final destination..."
            mv "$EXPECTED_SRA" "${OUT_DIR}/"

            echo "[INFO] Cleaning up temporary run workspace..."
            rm -rf "${TMP_DOWNLOAD_DIR}/${sample_id:?}"
            echo "[SUCCESS] Finished processing $sample_id"
        else
            echo "[ERROR] Prefetch claimed success, but file was not found at: $EXPECTED_SRA"
        fi
    else
        echo "[ERROR] Prefetch failed to download Run ID: $sample_id"
    fi

done

# 6. Final Cleanup -> (clean the directories where the files were downloaded)
if [ -d "$TMP_DOWNLOAD_DIR" ]; then
    rmdir "$TMP_DOWNLOAD_DIR" 2>/dev/null || true
fi

echo "=========================================================="
echo "[INFO] Processing complete! Files are saved in: $OUT_DIR"
echo "=========================================================="