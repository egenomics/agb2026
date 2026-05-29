#!/usr/bin/env python3
import sys
# For this, we need to load the container in the nextflow script
import pandas as pd

# Ensure proper arguments are passed
if len(sys.argv) != 3:
    print("Usage: clean_multiqc.py <input_multiqc_fastqc.txt> <output_file.tsv>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Load the MultiQC FastQC summary file
df = pd.read_csv(input_file, sep='\t')
out_df = pd.DataFrame()

# Retrieve the sra_id
out_df['sample-id'] = df['Sample']

# Parameter: Length
out_df['length'] = df['avg_sequence_length']
def eval_length(l):
    if 120 <= l <= 320: return 'pass'
    elif l < 120: return 'few'
    else: return 'more'
out_df['length_interp'] = out_df['length'].apply(eval_length)

# Parameter: Deduplication
out_df['deduplicated'] = df['total_deduplicated_percentage']
def eval_dedup(d):
    if 80 <= d <= 95: return 'pass'
    elif d < 80: return 'few'
    else: return 'more'
out_df['deduplicated_interp'] = out_df['deduplicated'].apply(eval_dedup)

# Parameter: %GC
out_df['%GC'] = df['%GC']
def eval_gc(gc):
    if 40 <= gc <= 60: return 'pass'
    elif gc < 40: return 'few'
    else: return 'more'
out_df['gc_interp'] = out_df['%GC'].apply(eval_gc)

# Parameter: Quality score
out_df['quality_score_status'] = df['per_sequence_quality_scores']
def eval_quality(q):
    if str(q).lower() == 'pass': return 'pass'
    else: return 'fail'
out_df['quality_interp'] = out_df['quality_score_status'].apply(eval_quality)

# Save the results in the output file
out_df.to_csv(output_file, sep="\t", index=False)
