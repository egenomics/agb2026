#!/usr/bin/bash
#SBATCH --job-name=diversity_pipeline
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/nextflow_%j.log

set -euo pipefail
mkdir -p logs

module load conda
module load nextflow

nextflow run diversity_analysis.nf \
  --data_dir ./data \
  --outdir ./results
