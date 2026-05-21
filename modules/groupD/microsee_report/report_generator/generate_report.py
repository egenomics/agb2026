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
import datetime
import hashlib
import logging
import platform
import subprocess
import sys
import time
from pathlib import Path

if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# Heavy imports (numpy, pandas) are deferred to main() so that
# --help and argument errors print instantly without a cold-start delay.

log = logging.getLogger(__name__)

# ── Version ───────────────────────────────────────────────────────────────────

_VERSION: str | None = None


def _get_version() -> str:
    global _VERSION
    if _VERSION is None:
        try:
            from importlib.metadata import version

            _VERSION = version("microsee-report")
        except Exception:
            _VERSION = "unknown"
    return _VERSION


# ── Logging ───────────────────────────────────────────────────────────────────


def _configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[MicroSee] %(levelname)s %(message)s",
        stream=sys.stderr,
    )


# ── CLI ───────────────────────────────────────────────────────────────────────

_EPILOG = """\
examples:
  microsee-report \\
      --feature-table feature-table.tsv \\
      --taxonomy      taxonomy.tsv      \\
      --metadata      metadata.tsv      \\
      --alpha         alpha-diversity.tsv

  # Cohort report + one HTML per patient:
  microsee-report ... --mode all --output results/report.html

  # Demo with bundled fixture data (from repo root):
  microsee-report \\
      --feature-table modules/groupD/microsee_report/tests/data/feature-table.tsv \\
      --taxonomy      modules/groupD/microsee_report/tests/data/taxonomy.tsv      \\
      --metadata      modules/groupD/microsee_report/tests/data/metadata.tsv      \\
      --alpha         modules/groupD/microsee_report/tests/data/alpha-diversity.tsv
"""


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "MicroSee — generate a self-contained offline HTML microbiome report.\n"
            "Reads QIIME2 TSV exports and writes one HTML file (~5 MB) with 30+ "
            "interactive Plotly charts. No internet required — Plotly.js is embedded."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_EPILOG,
    )
    p.add_argument("--version", action="version", version=f"microsee-report {_get_version()}")
    p.add_argument(
        "--feature-table",
        required=True,
        metavar="TSV",
        help=(
            "QIIME2 feature-table.tsv export. "
            "Columns = sample IDs; rows = ASV/OTU IDs; values = read counts."
        ),
    )
    p.add_argument(
        "--taxonomy",
        required=True,
        metavar="TSV",
        help=(
            "QIIME2 taxonomy.tsv export (Feature ID | Taxon | Confidence). "
            "Taxon strings must use SILVA or Greengenes format with 'f__' family prefixes."
        ),
    )
    p.add_argument(
        "--metadata",
        required=True,
        metavar="TSV",
        help=(
            "Sample metadata TSV. Required column: 'sample-id'. "
            "Recommended: 'subject' (patient ID), 'group' (treatment arm), 'timepoint' (T0/T84). "
            "Optional clinical columns: 'sixmwt' (metres), 'il18' (pg/mL)."
        ),
    )
    p.add_argument(
        "--alpha",
        metavar="TSV",
        default=None,
        help=(
            "Merged alpha-diversity TSV with per-sample metrics: "
            "shannon_entropy, simpson, faith_pd, pielou_evenness, observed_features. "
            "Strongly recommended — without this, diversity metrics are estimated from "
            "family-level abundances, which underestimates richness."
        ),
    )
    p.add_argument(
        "--distance-matrix",
        metavar="TSV",
        default=None,
        help=(
            "QIIME2 distance-matrix.tsv (square symmetric matrix, e.g. Bray-Curtis). "
            "When provided, this is used for PCoA / NMDS / PERMANOVA instead of the "
            "internally computed Bray-Curtis matrix. Optional."
        ),
    )
    p.add_argument(
        "--output",
        "-o",
        metavar="HTML",
        default="microsee_report.html",
        help=(
            "Output HTML file path (default: microsee_report.html). "
            "Parent directories are created automatically. "
            "In --mode all, per-patient files are written alongside this file."
        ),
    )
    p.add_argument(
        "--mode",
        choices=["cohort", "patient", "all"],
        default="cohort",
        help=(
            "cohort  — one combined report for all samples (default); "
            "patient — one HTML per patient (no cohort report); "
            "all     — cohort report + one HTML per patient."
        ),
    )
    return p.parse_args()


# ── Provenance ────────────────────────────────────────────────────────────────


class InputError(Exception):
    """Raised by _read_file when a CLI input cannot be read."""


def _read_file(path: str, label: str) -> str:
    p = Path(path)
    if not p.exists():
        raise InputError(
            f"{label} file not found: {path}\n  Tip: check the path is correct and the file exists."
        )
    try:
        # utf-8-sig strips the UTF-8 BOM that Excel and some Windows tools prepend.
        # Plain UTF-8 files are read identically — no downside.
        return p.read_text(encoding="utf-8-sig")
    except OSError as exc:
        raise InputError(f"Cannot read {label} file {path}: {exc}") from exc


def _file_sha256(path: str) -> str:
    """Return hex SHA-256 of a file. Used for provenance traceability."""
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return "unavailable"


def _get_git_hash() -> str:
    """Return short git commit hash, or 'unknown' if git is unavailable."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass
    return "unknown"


def _get_dep_versions() -> dict[str, str]:
    """Return installed versions of runtime dependencies."""
    from importlib.metadata import version as _v

    deps = {}
    for pkg in ("pandas", "numpy", "pydantic"):
        try:
            deps[pkg] = _v(pkg)
        except Exception:
            deps[pkg] = "unknown"
    return deps


def _build_provenance(args: argparse.Namespace) -> dict:
    """Collect full reproducibility metadata for this report run."""
    now_utc = datetime.datetime.now(datetime.UTC)
    input_files: dict[str, dict[str, str] | None] = {
        "feature_table": {
            "path": str(Path(args.feature_table).resolve()),
            "sha256": _file_sha256(args.feature_table),
        },
        "taxonomy": {
            "path": str(Path(args.taxonomy).resolve()),
            "sha256": _file_sha256(args.taxonomy),
        },
        "metadata": {
            "path": str(Path(args.metadata).resolve()),
            "sha256": _file_sha256(args.metadata),
        },
        "alpha": (
            {
                "path": str(Path(args.alpha).resolve()),
                "sha256": _file_sha256(args.alpha),
            }
            if args.alpha
            else None
        ),
        "distance_matrix": (
            {
                "path": str(Path(args.distance_matrix).resolve()),
                "sha256": _file_sha256(args.distance_matrix),
            }
            if args.distance_matrix
            else None
        ),
    }
    return {
        "generated_at": now_utc.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "software": "microsee-report",
        "version": _get_version(),
        "git_commit": _get_git_hash(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "dependencies": _get_dep_versions(),
        "input_files": input_files,
    }


# ── Main ──────────────────────────────────────────────────────────────────────


def main() -> None:
    _configure_logging()
    args = parse_args()
    _t_start = time.perf_counter()

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
        log.error("%s", exc)
        sys.exit(1)
    except ValueError as exc:
        log.error("Input validation failed: %s", exc)
        log.error("Check that all TSV files are valid QIIME2 exports with matching sample IDs.")
        sys.exit(1)

    for w in result.warnings:
        log.warning("%s", w)

    # ── Performance safeguards ────────────────────────────────────────────────
    if result.n_samples > 500:
        log.warning(
            "Large cohort: %d samples detected. "
            "Chart computation may take several minutes and the output HTML may exceed 20 MB. "
            "Consider using --mode cohort and skipping per-patient reports.",
            result.n_samples,
        )
    elif result.n_samples > 200:
        log.warning(
            "Large cohort: %d samples detected. Output HTML may exceed 10 MB.",
            result.n_samples,
        )

    if result.n_taxa > 100:
        log.warning(
            "%d taxa families detected. Distance matrix and heatmap operations will be slow. "
            "Consider filtering low-abundance taxa before generating the report.",
            result.n_taxa,
        )

    # ── Provenance ────────────────────────────────────────────────────────────
    provenance = _build_provenance(args)
    log.info(
        "Provenance: version=%s  git=%s  generated=%s",
        provenance["version"],
        provenance["git_commit"],
        provenance["generated_at"],
    )

    out = Path(args.output)
    # Ensure the output directory exists — useful when --output is a deep path.
    out.parent.mkdir(parents=True, exist_ok=True)

    if args.mode in ("cohort", "all"):
        log.info("Computing chart data...")
        _t_charts = time.perf_counter()
        chart_data = compute_chart_data(result, distance_matrix=distance_matrix)
        log.info("Chart data computed in %.1f s", time.perf_counter() - _t_charts)
        chart_data["meta"]["mode"] = args.mode
        chart_data["meta"]["provenance"] = provenance
        if args.mode == "all":
            stem = out.stem
            suffix = out.suffix or ".html"
            chart_data["meta"]["patient_nav"] = [
                {
                    "label": pid,
                    "href": f"{stem}_{pid.replace('/', '_').replace(' ', '_')}{suffix}",
                }
                for pid in sorted({r.patient for r in result.rows})
            ]

        log.info("Rendering HTML...")
        try:
            html = render_html(chart_data)
        except FileNotFoundError as exc:
            log.error("%s", exc)
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
                    provenance=provenance,
                )
            except FileNotFoundError as exc:
                log.error("%s", exc)
                sys.exit(1)
            safe_pid = patient_id.replace("/", "_").replace(" ", "_")
            patient_out = parent / f"{stem}_{safe_pid}{suffix}"
            patient_out.write_text(patient_html, encoding="utf-8")
            log.info("  Patient report → %s", patient_out.name)

    log.info("Total pipeline time: %.1f s", time.perf_counter() - _t_start)


if __name__ == "__main__":
    main()
