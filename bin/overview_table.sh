#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(dirname "$(readlink -f "$0")") 

METADATA=$1
PCA=$2
ZSCORES=$3
GENES_COUNTS=$4
OUTDIR=$5

Rscript ${SCRIPT_DIR}/overview_table.R \
  --metadata $METADATA \
  --pca $PCA \
  --zscored $ZSCORES \
  --counts $GENES_COUNTS \
  --outdir $OUTDIR