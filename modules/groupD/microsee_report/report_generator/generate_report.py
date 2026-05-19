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
        [--distance-matrix distance-matrix.tsv] \\
        [--output       microsee_report.html]

Usage (direct):
    python generate_report.py --feature-table ...
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from report_generator.charts import compute_chart_data, render_html, render_patient_html
from report_generator.charts.individual import build_patient_radar_profiles
from report_generator.models import DistanceMatrixResult
from report_generator.parsers import (
    integrate,
    parse_alpha_diversity,
    parse_distance_matrix,
    parse_feature_table,
    parse_metadata,
    parse_taxonomy,
)

log = logging.getLogger(__name__)


def _configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[MicroSee] %(levelname)s %(message)s",
        stream=sys.stderr,
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate MicroSee self-contained HTML report")
    p.add_argument(
        "--feature-table", required=True, metavar="TSV", help="QIIME2 feature-table.tsv export"
    )
    p.add_argument("--taxonomy", required=True, metavar="TSV", help="QIIME2 taxonomy.tsv export")
    p.add_argument("--metadata", required=True, metavar="TSV", help="QIIME2 metadata.tsv export")
    p.add_argument(
        "--alpha",
        metavar="TSV",
        default=None,
        help="QIIME2 alpha-diversity.tsv export (strongly recommended)",
    )
    p.add_argument(
        "--distance-matrix",
        metavar="TSV",
        default=None,
        help="QIIME2 distance-matrix.tsv (e.g. Bray-Curtis); optional",
    )
    p.add_argument(
        "--output",
        "-o",
        metavar="HTML",
        default="microsee_report.html",
        help="Output HTML file (default: microsee_report.html)",
    )
    p.add_argument(
        "--mode",
        choices=["cohort", "patient", "all"],
        default="cohort",
        help=(
            "cohort  — one combined report (default); "
            "patient — one HTML per patient; "
            "all     — cohort report + one per patient"
        ),
    )
    return p.parse_args()


class InputError(Exception):
    """Raised by _read_file when a CLI input cannot be read."""


def _read_file(path: str, label: str) -> str:
    p = Path(path)
    if not p.exists():
        raise InputError(f"{label} file not found: {path}")
    try:
        return p.read_text(encoding="utf-8")
    except OSError as exc:
        raise InputError(f"Cannot read {label} file {path}: {exc}") from exc


def main() -> None:
    _configure_logging()
    args = parse_args()

    distance_matrix: DistanceMatrixResult | None = None

    try:
        log.info("Parsing feature table:   %s", args.feature_table)
        feat = parse_feature_table(_read_file(args.feature_table, "feature-table"))

        log.info("Parsing taxonomy:         %s", args.taxonomy)
        tax = parse_taxonomy(_read_file(args.taxonomy, "taxonomy"))

        log.info("Parsing metadata:         %s", args.metadata)
        meta = parse_metadata(_read_file(args.metadata, "metadata"))

        alpha = None
        if args.alpha:
            log.info("Parsing alpha diversity:  %s", args.alpha)
            alpha = parse_alpha_diversity(_read_file(args.alpha, "alpha-diversity"))

        if args.distance_matrix:
            log.info("Parsing distance matrix:  %s", args.distance_matrix)
            distance_matrix = parse_distance_matrix(
                _read_file(args.distance_matrix, "distance-matrix"),
            )

        log.info("Integrating data...")
        result = integrate(feat, tax, meta, alpha)

    except InputError as exc:
        log.error("[MicroSee] %s", exc)
        sys.exit(1)
    except ValueError as exc:
        log.error("[MicroSee] Input validation failed: %s", exc)
        log.error("Check that all TSV files are valid QIIME2 exports with matching sample IDs.")
        sys.exit(1)

    for w in result.warnings:
        log.warning("%s", w)

    out = Path(args.output)

    if args.mode in ("cohort", "all"):
        log.info("Computing chart data...")
        chart_data = compute_chart_data(result, distance_matrix=distance_matrix)

        log.info("Rendering HTML...")
        try:
            html = render_html(chart_data)
        except FileNotFoundError as exc:
            log.error("[MicroSee] %s", exc)
            sys.exit(1)
        out.write_text(html, encoding="utf-8")
        log.info("Report written → %s  (%.1f KB)", out.resolve(), out.stat().st_size / 1024)

    if args.mode in ("patient", "all"):
        patients = sorted({r.patient for r in result.rows})
        stem = out.stem
        suffix = out.suffix or ".html"
        parent = out.parent
        log.info("Generating per-patient reports for %d patients ...", len(patients))
        rows_dump = [r.model_dump() for r in result.rows]
        radar_profiles = build_patient_radar_profiles(rows_dump, result.taxa)
        for patient_id in patients:
            try:
                patient_html = render_patient_html(
                    patient_id,
                    result,
                    radar_profiles=radar_profiles,
                )
            except FileNotFoundError as exc:
                log.error("[MicroSee] %s", exc)
                sys.exit(1)
            safe_pid = patient_id.replace("/", "_").replace(" ", "_")
            patient_out = parent / f"{stem}_{safe_pid}{suffix}"
            patient_out.write_text(patient_html, encoding="utf-8")
            log.info("  Patient report → %s", patient_out.name)


if __name__ == "__main__":
    main()
