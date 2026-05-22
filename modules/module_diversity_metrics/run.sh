#!/usr/bin/bash
#SBATCH --job-name=diversity_pipeline
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/nextflow_%j.log

set -euo pipefail

mkdir -p logs   # make sure logs dir exists before SLURM tries to write to it

module load conda
module load nextflow 

nextflow run diversity_analysis.nf \
  --data_dir ./data \
  --outdir results \
  --sampling_depth 1103
