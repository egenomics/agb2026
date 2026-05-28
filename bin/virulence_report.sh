#!/bin/bash


set -euo pipefail

ASV_TABLE_PATH=$1
TAXONOMY_PATH=$2
METADATA_PATH=$3

python3 virulence_report.py "$ASV_TABLE_PATH" "$TAXONOMY_PATH" "$METADATA_PATH"