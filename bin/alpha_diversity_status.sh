#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(dirname "$(readlink -f "$0")") 

METADATA=$1
FAITH=$2
OBSERVED=$3
SHANNON=$4
SIMPSON=$5
OUTDIR=$6

Rscript ${SCRIPT_DIR}/alpha_diversity_status.R \
  --metadata $METADATA \
  --faith $FAITH \
  --observed $OBSERVED \
  --shannon $SHANNON \
  --simpson $SIMPSON \
  --outdir $OUTDIR