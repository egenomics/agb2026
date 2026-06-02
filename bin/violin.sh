#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

ASV_TABLE_PATH=$1
TAXONOMY_PATH=$2
METADATA_PATH=$3
OUTDIR=$4

python3 ${SCRIPT_DIR}/Violin_plot.py --asv-table "$ASV_TABLE_PATH" --taxonomy "$TAXONOMY_PATH" --metadata "$METADATA_PATH" --outdir $OUTDIR