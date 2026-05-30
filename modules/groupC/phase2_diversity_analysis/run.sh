#!/usr/bin/bash
 
set -euo pipefail
mkdir -p logs
 
nextflow run main.nf \
  --data_dir modules/module_diversity_metrics/data \
  --outdir   results \
  -with-report logs/report.html \
  -with-trace  logs/trace.txt