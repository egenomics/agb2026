#!/usr/bin/bash
#SBATCH --job-name=contamination_filter_pipeline
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/nextflow_%j.log

set -euo pipefail

mkdir -p logs

module load conda
module load nextflow

nextflow run contamination_filter.nf \
  --data_dir ./data \
  --outdir   results