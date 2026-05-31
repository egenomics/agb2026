#!/bin/bash


set -euo pipefail

ASV_TABLE_PATH=$1
TAXONOMY_PATH=$2
METADATA_PATH=$3

python3 virulence_report.py \
    --asv-table "$ASV_TABLE_PATH" \
    --taxonomy "$TAXONOMY_PATH" \
    --metadata "$METADATA_PATH" \
    --outdir "Virulence_analysis"