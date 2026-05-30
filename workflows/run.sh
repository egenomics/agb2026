#!/usr/bin/bash
#SBATCH --job-name=diversity_pipeline
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/nextflow_%j.log

set -euo pipefail

if [[ $(hostname) == *"login"* ]]; then
    echo "ERROR: Do not run this script on a login node. Use: sbatch run.sh"
    exit 1
fi

mkdir -p logs

module load apptainer
module load nextflow

if ! command -v apptainer &> /dev/null; then
    echo "ERROR: apptainer not found after module load"
    exit 1
fi

apptainer --version

nextflow run groupC.nf \
  --data_dir ../modules/groupC/phase2_diversity_analysis/modules/module_diversity_metrics/data \
  --outdir   results \
  -work-dir  /data/upfagb/u269208/agb2026/workflows/work \
  -with-report logs/report.html \
  -with-trace  logs/trace.txt \
  -resume
