#!/usr/bin/env bash
set -euo pipefail

# This script must be run from the stress_tests directory:
# modules/groupC/phase1_system_validation/stress_tests

GT="../dataset_validation/ground_truth_PRJEB10949.tsv"
VM="../validation_metrics.py"
TAX="ASV_taxonomy.tsv"
OUTDIR="stress_results"

SCENARIOS=(
  ST00_baseline
  ST01_zero_count_sample
  ST02_low_depth
  ST04_single_taxon
  ST06b_contamination_biological
)

# Basic checks before running
[[ -f "$GT" ]] || { echo "ERROR: Ground truth not found: $GT"; exit 1; }
[[ -f "$VM" ]] || { echo "ERROR: validation_metrics.py not found: $VM"; exit 1; }
[[ -f "$TAX" ]] || { echo "ERROR: Taxonomy file not found: $TAX"; exit 1; }

mkdir -p "$OUTDIR"

for SC in "${SCENARIOS[@]}"; do
  ASV_TABLE="stress_inputs/$SC/asv_table.tsv"
  [[ -f "$ASV_TABLE" ]] || { echo "ERROR: ASV table not found: $ASV_TABLE"; exit 1; }

  echo "Running validation metrics for $SC"
  mkdir -p "$OUTDIR/$SC"

  python "$VM" "$GT" "$ASV_TABLE" "$TAX"

  # Move outputs if validation_metrics.py writes them in the current directory
  [[ -f detection_metrics.tsv ]] && mv detection_metrics.tsv "$OUTDIR/$SC/"
  [[ -f abundance_metrics.tsv ]] && mv abundance_metrics.tsv "$OUTDIR/$SC/"
done

echo "All validation metrics completed. Results saved in: $OUTDIR"
