#!/usr/bin/bash
#SBATCH --job-name=kraken2_pipeline
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/nextflow_%j.log

set -euo pipefail

mkdir -p logs

module load conda
module load nextflow

nextflow run kraken2.nf \
  --data_dir ./data \
  --outdir   results \
  --kraken2_db /data/upfagb/u269238/kraken2_db