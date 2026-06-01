#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

ASV_TABLE_PATH=$1
TAXONOMY_PATH=$2
METADATA_PATH=$3

python3 ${SCRIPT_DIR}/Violin_plot.py "$ASV_TABLE_PATH" "$TAXONOMY_PATH" "$METADATA_PATH"