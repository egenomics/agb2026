#!/usr/bin/env python3
"""
generate_report.py — MicroSee self-contained HTML report generator (CLI entry point).

Reads QIIME2 TSV exports, parses them, pre-computes all chart data, and writes
one HTML file with embedded Plotly.js charts.  No server needed — open in any browser.

Usage (after pip install -e .):
    microsee-report \\
        --feature-table feature-table.tsv \\
        --taxonomy      taxonomy.tsv      \\
        --metadata      metadata.tsv      \\
        [--alpha        alpha-diversity.tsv] \\
        [--output       microsee_report.html]

Usage (direct):
    python generate_report.py --feature-table ...

Chart logic lives in charts/:
    config.py       — colour constants, Plotly layout/config defaults
    utils.py        — shared colour helpers (group/taxon palette)
    distances.py    — Bray-Curtis, Jaccard, PCoA, average-linkage clustering
    taxonomy.py     — stacked bar (27 filter variants), top-N, donut, sunburst
    alpha.py        — strip/box/violin, rarefaction, multi-metric, significance brackets
    beta.py         — PCoA (Bray-Curtis + Jaccard), NMDS, dendrogram, Δ heatmap
    individual.py   — slopegraph, stability bar, rank plot, radar, faceted small-multiples, trajectories
    comparative.py  — LFC bar, volcano (BH-FDR), ANCOM-style CLR, heatmap, correlation matrix
    clinical.py     — clinical slopegraph, Shannon correlation, taxa×clinical Spearman heatmap
    stats.py        — Wilcoxon/Mann-Whitney, longitudinal, LME trajectory, diversity table, PERMANOVA
    renderer.py     — compute_chart_data(), render_html()
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Allow running as a plain script (python generate_report.py) without installation.
# When installed via pip, this block is skipped because the package is on sys.path.
if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from report_generator.parsers import (
    parse_alpha_diversity,
    parse_feature_table,
    parse_metadata,
    parse_taxonomy,
    integrate,
)
from report_generator.charts import compute_chart_data, render_html
from report_generator.charts.renderer import render_patient_html

log = logging.getLogger(__name__)


def _configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[MicroSee] %(levelname)s %(message)s",
        stream=sys.stderr,
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate MicroSee self-contained HTML report")
    p.add_argument("--feature-table", required=True, metavar="TSV",
                   help="QIIME2 feature-table.tsv export")
    p.add_argument("--taxonomy",      required=True, metavar="TSV",
                   help="QIIME2 taxonomy.tsv export")
    p.add_argument("--metadata",      required=True, metavar="TSV",
                   help="QIIME2 metadata.tsv export")
    p.add_argument("--alpha",         metavar="TSV", default=None,
                   help="QIIME2 alpha-diversity.tsv export (optional)")
    p.add_argument("--output", "-o",  metavar="HTML",
                   default="microsee_report.html",
                   help="Output HTML file (default: microsee_report.html)")
    p.add_argument("--mode", choices=["cohort", "patient", "all"],
                   default="cohort",
                   help=(
                       "cohort  — one combined report (default); "
                       "patient — one HTML per patient; "
                       "all     — cohort report + one per patient"
                   ))
    return p.parse_args()


def main() -> None:
    _configure_logging()
    args = parse_args()

    def read(path: str) -> str:
        return Path(path).read_text(encoding="utf-8")

    log.info("Parsing feature table:   %s", args.feature_table)
    feat = parse_feature_table(read(args.feature_table))

    log.info("Parsing taxonomy:         %s", args.taxonomy)
    tax  = parse_taxonomy(read(args.taxonomy))

    log.info("Parsing metadata:         %s", args.metadata)
    meta = parse_metadata(read(args.metadata))

    alpha = None
    if args.alpha:
        log.info("Parsing alpha diversity:  %s", args.alpha)
        alpha = parse_alpha_diversity(read(args.alpha))

    log.info("Integrating data...")
    result = integrate(feat, tax, meta, alpha)

    for w in result.warnings:
        log.warning("%s", w)

    out = Path(args.output)

    if args.mode in ("cohort", "all"):
        log.info("Computing chart data...")
        chart_data = compute_chart_data(result)

        log.info("Rendering HTML...")
        html = render_html(chart_data)
        out.write_text(html, encoding="utf-8")
        log.info("Report written → %s  (%.1f KB)", out.resolve(), out.stat().st_size / 1024)

    if args.mode in ("patient", "all"):
        patients = sorted(set(r.patient for r in result.rows))
        stem     = out.stem if args.mode == "all" else out.stem
        suffix   = out.suffix or ".html"
        parent   = out.parent
        log.info("Generating per-patient reports for %d patients ...", len(patients))
        for patient_id in patients:
            patient_html = render_patient_html(patient_id, result)
            safe_pid     = patient_id.replace("/", "_").replace(" ", "_")
            patient_out  = parent / f"{stem}_{safe_pid}{suffix}"
            patient_out.write_text(patient_html, encoding="utf-8")
            log.info("  Patient report → %s", patient_out.name)


if __name__ == "__main__":
    main()
