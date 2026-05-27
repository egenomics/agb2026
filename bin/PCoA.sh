#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(dirname "$(readlink -f "$0")") 

METADATA=$1
BRAY=$2
W_UNIFRAC=$3
OUTDIR=$4

Rscript PCoA.R