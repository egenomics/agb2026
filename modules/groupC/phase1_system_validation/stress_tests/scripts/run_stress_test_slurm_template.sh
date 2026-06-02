#!/bin/bash
#SBATCH --job-name=agb_groupC_stress
#SBATCH --output=/data/upfagb/<USER>/agb_stress_tests/logs/%x_%j.out
#SBATCH --error=/data/upfagb/<USER>/agb_stress_tests/logs/%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

# Edit these paths before running.
REPO=/data/upfagb/<USER>/agb2026
OUTDIR=/data/upfagb/<USER>/agb_stress_tests/results
WORKDIR=/scratch/upfagb/<USER>/agb_stress_tests/work

mkdir -p "$OUTDIR" "$WORKDIR"

cd "$REPO"

echo "Started at: $(date)"
echo "Node: $(hostname)"
echo "Working directory: $PWD"

# Load modules if required by the cluster.
# module load nextflow
# module load singularity

# Example only. Replace this with the real Group C / full pipeline command
# once the interface is agreed with Group B.
#
# nextflow run main.nf \
#   -profile singularity,slurm \
#   --input modules/groupC/phase1_system_validation/stress_tests/stress_inputs/ST00_valid_subset/asv_table_subset.tsv \
#   --outdir "$OUTDIR/ST00_valid_subset" \
#   -work-dir "$WORKDIR/ST00_valid_subset"

echo "Template finished at: $(date)"
