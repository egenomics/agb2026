#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# build_report.sh — Group D MicroSee report wrapper
# AGB 2026
#
# Matches the Group D bin/ convention (see alpha_diversity_status.sh, PCoA.sh):
# the .nf process calls this; it stages the incoming plot directories into one
# input tree, then runs the generator, which scans recursively.
#
# Usage:
#   build_report.sh <outdir> <plot_dir> [<plot_dir> ...] [-- PATIENT_ID]
#
#   <outdir>      output directory for the report
#   <plot_dir>    one or more directories emitted by upstream plotting
#                 processes (exploratory/, explanatory_report/, volcano_plot/,
#                 beta_diversity_tree/, demographic_table/,
#                 clinical_association_map/, ...). Any number, any nesting.
#   PATIENT_ID    optional, after a literal `--`: build only this patient.
#
# The generator (build_album.py) recurses through whatever is staged:
#   - PNGs whose name contains an ERR id  -> per-patient report
#   - all other PNGs                      -> cohort (exploratory) report
#   - non-PNG files (.rds, .tsv, .qza)    -> ignored
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

OUTDIR=$1
shift

PATIENT_ID=""
PLOT_DIRS=()
while [ $# -gt 0 ]; do
    if [ "$1" = "--" ]; then
        shift
        PATIENT_ID=${1:-}
        break
    fi
    PLOT_DIRS+=("$1")
    shift
done

mkdir -p build_input
for d in "${PLOT_DIRS[@]}"; do
    [ -e "$d" ] || continue
    base=$(basename "$d")

    clean_base=$(echo "$base" | sed -E 's/\.[0-9]+$//')

    if [ -d "$d" ]; then
        mkdir -p "build_input/$clean_base"
        if [ "$(ls -A "$d")" ]; then
            cp -rL "$d"/. "build_input/$clean_base/" 2>/dev/null || true
        fi
    else
        # If it's a standalone file, copy it directly under the clean name
        cp -rL "$d" "build_input/$clean_base" 2>/dev/null || true
    fi
done

if [ -n "$PATIENT_ID" ]; then
    python3 "${SCRIPT_DIR}/build_album.py" \
        --input build_input --patient "$PATIENT_ID" --outdir "$OUTDIR"
else
    python3 "${SCRIPT_DIR}/build_album.py" \
        --input build_input --outdir "$OUTDIR"
fi
