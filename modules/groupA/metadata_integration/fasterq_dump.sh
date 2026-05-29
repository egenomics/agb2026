#!/bin/bash

# Prepare a conda enviorment to run the function (or your own alternative)
# conda create -n fasterq_dump_env -y 
conda activate fasterq_dump_env
# conda install bioconda::sra-tools -y 

# 1. Help function
show_help() {
    cat << EOF
Usage: ${0##*/} [-h] -i INPUT_DIR -o OUTPUT_DIR -s SAMPLE_SHEET

Description: Automates fasterq-dump for a list of .sra files defined in a metadata TSV.
             Dynamically locates 'sample-id', reads data via Conda safely, and 
             automatically handles SE/PE layouts using --split-3.

OPTIONS:
    -h         Show this help message.
    -i DIR     Path to the input directory containing the raw .sra files.
    -o DIR     Path to the base output directory for fastq files.
    -s FILE    Path to the metadata sample sheet file (TSV format).
EOF
}

# 2. Parse Arguments
IN_DIR=""
OUT_DIR=""
SAMPLE_SHEET=""

while getopts "hi:o:s:" opt; do
    case "$opt" in
        h) show_help; exit 0 ;;
        i) IN_DIR=$OPTARG ;;
        o) OUT_DIR=$OPTARG ;;
        s) SAMPLE_SHEET=$OPTARG ;;
        *) show_help; exit 1 ;;
    esac
done

# Ensure mandatory arguments are provided
if [ -z "$IN_DIR" ] || [ -z "$OUT_DIR" ] || [ -z "$SAMPLE_SHEET" ]; then
    echo "Error: Missing required arguments (-i, -o, and -s are mandatory)."
    show_help
    exit 1
fi

# Ensure input directory and sample sheet exist
if [ ! -d "$IN_DIR" ]; then
    echo "Error: Input directory '$IN_DIR' not found."
    exit 1
fi

if [ ! -f "$SAMPLE_SHEET" ]; then
    echo "Error: Sample sheet file '$SAMPLE_SHEET' not found."
    exit 1
fi

# 3. Setup Outputs & Conda Environment
FINAL_OUT_DIR="${OUT_DIR}/inter_data/seqs/splitted"
mkdir -p "$FINAL_OUT_DIR"

# 4. Find the 'sample-id' column
DELIM=$'\t'
header=$(head -n 1 "$SAMPLE_SHEET" | tr -d '\r' | tr -d '"')
col_idx=$(echo "$header" | tr "$DELIM" '\n' | grep -nx "sample-id" | cut -d: -f1)

if [ -z "$col_idx" ]; then
    echo "Error: Could not find 'sample-id' column in the metadata header."
    exit 1
fi

echo "[INFO] Successfully mapped 'sample-id' to column index: $col_idx"

# 5. Iterate over the TSV and extract
echo "[INFO] Starting FASTQ extraction..."

# tail -n +2 skips the header row
tail -n +2 "$SAMPLE_SHEET" | cut -d"$DELIM" -f"$col_idx" | while read -r sample_raw; do
    
    # Clean up hidden carriage returns or wrapping quotes
    sample_id=$(echo "$sample_raw" | tr -d '\r' | tr -d '"' | xargs)
    
    # Skip empty lines
    if [ -z "$sample_id" ]; then
        continue
    fi

    SRA_FILE="${IN_DIR}/${sample_id}.sra"

    echo "=========================================================="
    echo "[INFO] Extracting Sequence for Run: $sample_id"
    echo "=========================================================="

    # Check if the raw SRA file actually exists before running fasterq-dump
    if [ ! -f "$SRA_FILE" ]; then
        echo "[WARNING] Target file not found: $SRA_FILE. Skipping."
        continue
    fi

    # Run fasterq-dump using 'conda run' to avoid shell activation bugs.
    # --split-3 automatically detects SE or PE structures and writes files accordingly.
    # </dev/null keeps the stream clean for our while-read loop.
    if fasterq-dump "$SRA_FILE" \
        -O "$FINAL_OUT_DIR" \
        --split-3 </dev/null; then
        echo "[SUCCESS] Finished extracting $sample_id"
    else
        echo "[ERROR] fasterq-dump execution failed for $sample_id"
    fi

done

echo "=========================================================="
echo "[INFO] Extraction complete! Fastq files saved in: $FINAL_OUT_DIR"
echo "=========================================================="