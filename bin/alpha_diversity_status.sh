#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(dirname "$(readlink -f "$0")") 

METADATA=$1
BRAY=$2
UNIFRAC=$3
FAITH=$4
OBSERVED=$5
SHANNON=$6
SIMPSON=$7
OUTDIR=$8

Rscript ${SCRIPT_DIR}/alpha_diversity_status.R \
  --metadata $1 \
  --bray $2 \
  --unifrac $3 \
  --faith $4 \
  --observed $5 \
  --shannon $6 \
  --simpson $7 \
  --outdir $8