#!/usr/bin/env python
import sys
# For this, we need to load the container in the nextflow script
import pandas as pd

# Basic error handling to ensure Nextflow passes the right number of arguments
if len(sys.argv) != 4:
    print("Usage: generate_sample_sheet_B.py <id_list.txt> <trimmed_dir_path> <output_file.csv>")
    sys.exit(1)

# Map the arguments passed by Nextflow to variables
id_list_file = sys.argv[1]
trimmed_dir  = sys.argv[2]
output_csv   = sys.argv[3]

# Read the successful sample IDs
with open(id_list_file, 'r') as f:
    samples = [line.strip() for line in f if line.strip()]

# Build the data rows for each of the samples
data = {
    'sample_id': samples,
    'directory': [trimmed_dir] * len(samples)
}

# Create the DataFrame and save as a clean CSV file
df = pd.DataFrame(data)
df.to_csv(output_csv, index=False)