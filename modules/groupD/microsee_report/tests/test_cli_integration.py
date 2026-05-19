"""
Integration smoke tests — HTML output and CLI entry point.

Skipped by default (``pytest -m integration``). Uses a 4-sample inline dataset so
GitHub Actions finishes in ~1–3 minutes, not 20+.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from report_generator.charts import compute_chart_data, render_html, render_patient_html
from report_generator.parsers import (
    integrate,
    parse_alpha_diversity,
    parse_feature_table,
    parse_metadata,
    parse_taxonomy,
)

# Minimal inline fixtures (same as test_charts.py) — fast on CI runners.
FEATURE_TABLE = """\
#OTU ID\tS1\tS2\tS3\tS4
ASV1\t100\t200\t50\t80
ASV2\t300\t150\t400\t200
ASV3\t50\t50\t250\t100
"""
TAXONOMY = """\
Feature ID\tTaxon\tConfidence
ASV1\td__Bacteria;p__Firmicutes;f__Lachnospiraceae\t0.99
ASV2\td__Bacteria;p__Bacteroidota;f__Bacteroidaceae\t0.98
ASV3\td__Bacteria;p__Firmicutes;f__Ruminococcaceae\t0.97
"""
METADATA = """\
sample-id\tsubject\tgroup\ttimepoint
S1\tPat1\tEAA\tT0
S2\tPat1\tEAA\tT84
S3\tPat2\tEAA\tT0
S4\tPat2\tEAA\tT84
"""
ALPHA = """\
sample-id\tshannon_entropy\tsimpson\tfaith_pd
S1\t1.5\t0.75\t8.2
S2\t1.8\t0.82\t9.1
S3\t1.3\t0.70\t7.8
S4\t1.6\t0.78\t8.5
"""

PLOTLY_JS = Path(__file__).resolve().parents[1] / "report_generator" / "charts" / "plotly.min.js"
MIN_REPORT_BYTES = 1_000_000  # Plotly.js alone is ~4.3 MB
MIN_PATIENT_BYTES = 50_000


@pytest.fixture(scope="module")
def integrated_result():
    feat = parse_feature_table(FEATURE_TABLE)
    tax = parse_taxonomy(TAXONOMY)
    meta = parse_metadata(METADATA)
    alpha = parse_alpha_diversity(ALPHA)
    return integrate(feat, tax, meta, alpha)


@pytest.fixture(scope="module")
def cohort_html(integrated_result) -> str:
    """Full cohort HTML as a string — no disk I/O."""
    return render_html(compute_chart_data(integrated_result))


@pytest.mark.integration
def test_plotly_bundle_present() -> None:
    assert PLOTLY_JS.is_file(), (
        f"Missing bundled Plotly.js: {PLOTLY_JS}\n"
        "Commit report_generator/charts/plotly.min.js for offline HPC use."
    )
    assert PLOTLY_JS.stat().st_size > 1_000_000


@pytest.mark.integration
class TestCohortReport:
    def test_is_html(self, cohort_html: str) -> None:
        assert "<html" in cohort_html.lower()

    def test_size_above_minimum(self, cohort_html: str) -> None:
        size = len(cohort_html.encode())
        assert size >= MIN_REPORT_BYTES, (
            f"Cohort report too small: {size:,} bytes "
            f"(expected ≥ {MIN_REPORT_BYTES:,}). Plotly.js may not be embedded."
        )

    def test_plotly_charts_present(self, cohort_html: str) -> None:
        assert "Plotly.newPlot" in cohort_html

    def test_no_external_cdn(self, cohort_html: str) -> None:
        # plotly.min.js is embedded inline; its minified source may contain the
        # string "cdn.plot.ly" internally. The real check is: no external <script src=>.
        assert 'src="https://cdn.plot.ly' not in cohort_html


@pytest.mark.integration
class TestPatientReports:
    def test_patient_html_in_process(self, integrated_result) -> None:
        patient_id = integrated_result.rows[0].patient
        html = render_patient_html(patient_id, integrated_result)
        assert "<html" in html.lower()
        assert len(html) >= MIN_PATIENT_BYTES
        assert "Plotly.newPlot" in html


@pytest.mark.integration
class TestCliSubprocess:
    """CLI entry-point smoke — uses missing file so it exits before any heavy computation."""

    def test_missing_file_exits_nonzero(self, tmp_path: Path) -> None:
        result = subprocess.run(
            [
                sys.executable, "-m", "report_generator.generate_report",
                "--feature-table", str(tmp_path / "missing.tsv"),
                "--taxonomy", str(tmp_path / "taxonomy.tsv"),
                "--metadata", str(tmp_path / "metadata.tsv"),
                "--output", str(tmp_path / "report.html"),
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert result.returncode != 0
