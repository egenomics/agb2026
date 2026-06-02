#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(dirname "$(readlink -f "$0")") 

PCA=$1
METADATA=$2
OUTDIR=$3

Rscript ${SCRIPT_DIR}/PCoA.R --metadata $METADATA --pca $PCA --outdir $OUTDIR