#!/usr/bin/bash

set -euo pipefail
mkdir -p logs

nextflow run select_rarefaction.nf \
  --diversity_metrics_dir ../module_diversity_metrics \
  --outdir ./results \
  -with-report logs/report.html \
  -with-trace  logs/trace.txt
