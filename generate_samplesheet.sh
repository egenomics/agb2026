#!/bin/bash

# 1. Check if a directory path was provided
if [ -z "$1" ]; then
    echo "Usage: $0 /path/to/data"
    exit 1
fi

# 2. Get the absolute path so Nextflow doesn't get lost
DATA_DIR=$(realpath "$1")
OUTPUT_CSV="samplesheet.csv"

# 3. Write the exact header (4 columns — matches what test_groupB.nf parses)
echo "sample,fastq_1,fastq_2,strandedness" > "$OUTPUT_CSV"

# 4. Find all _trim.fastq.gz files
find "$DATA_DIR" -type f -name "*_trim.fastq.gz" | sort | while read -r FQ_FILE; do

    # 5. Extract the base sample name (e.g., ERR1074192)
    # %%_trim.* strips everything from "_trim." to the end — handles both
    # _trim.fastq and _trim.fastq.gz cleanly.
    BASENAME=$(basename "$FQ_FILE")
    SAMPLE_NAME="${BASENAME%%_trim.*}"

    # 6. Write to CSV (',,' leaves fastq_2 empty for single-end reads)
    echo "${SAMPLE_NAME},${FQ_FILE},,forward" >> "$OUTPUT_CSV"

done

N_SAMPLES=$(($(wc -l < "$OUTPUT_CSV") - 1))
echo "Success! Samplesheet generated at: $(pwd)/$OUTPUT_CSV  ($N_SAMPLES samples)"
