"""
test_cli_integration.py — end-to-end smoke tests for the microsee-report CLI.

Default ``pytest`` skips these (marker: integration). The session fixture runs
``--mode cohort`` only (~1–2 min). Per-patient HTML is checked in-process.
"""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

import pytest

from report_generator.charts import render_patient_html
from report_generator.parsers import (
    integrate,
    parse_alpha_diversity,
    parse_feature_table,
    parse_metadata,
    parse_taxonomy,
)

DATA = Path(__file__).parent / "data"
FEATURE_TABLE = DATA / "feature-table.tsv"
TAXONOMY      = DATA / "taxonomy.tsv"
METADATA      = DATA / "metadata.tsv"
ALPHA         = DATA / "alpha-diversity.tsv"

MIN_REPORT_BYTES  = 1_000_000  # 1 MB — Plotly.js alone is ~4.3 MB
MIN_PATIENT_BYTES = 50_000     # in-process patient report (no full cohort embed)


def _run(output: Path, *extra: str, timeout: int = 120) -> subprocess.CompletedProcess[str]:
    t0 = time.monotonic()
    result = subprocess.run(
        [
            sys.executable, "-m", "report_generator.generate_report",
            "--feature-table", str(FEATURE_TABLE),
            "--taxonomy",      str(TAXONOMY),
            "--metadata",      str(METADATA),
            "--alpha",         str(ALPHA),
            "--output",        str(output),
            *extra,
        ],
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    elapsed = time.monotonic() - t0
    if elapsed > 60:
        sys.stderr.write(f"[integration] CLI finished in {elapsed:.1f}s\n")
    return result


@pytest.fixture(scope="session")
def cohort_output(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """One subprocess: cohort report only (not --mode all — that writes 12+ huge HTML files)."""
    d = tmp_path_factory.mktemp("cohort")
    out = d / "microsee_report.html"
    result = _run(out, "--mode", "cohort", timeout=120)
    assert result.returncode == 0, f"CLI failed ({result.returncode}):\n{result.stderr}"
    return out


@pytest.fixture(scope="session")
def integrated_result():
    """Parsed fixture data for in-process patient render (no subprocess)."""
    feat = parse_feature_table(FEATURE_TABLE.read_text(encoding="utf-8"))
    tax = parse_taxonomy(TAXONOMY.read_text(encoding="utf-8"))
    meta = parse_metadata(METADATA.read_text(encoding="utf-8"))
    alpha = parse_alpha_diversity(ALPHA.read_text(encoding="utf-8"))
    return integrate(feat, tax, meta, alpha)


@pytest.mark.integration
class TestCohortReport:
    def test_file_exists(self, cohort_output: Path) -> None:
        assert cohort_output.exists()

    def test_is_html(self, cohort_output: Path) -> None:
        content = cohort_output.read_text(encoding="utf-8", errors="replace")
        assert "<html" in content.lower()

    def test_size_above_minimum(self, cohort_output: Path) -> None:
        size = cohort_output.stat().st_size
        assert size >= MIN_REPORT_BYTES, (
            f"Cohort report too small: {size:,} bytes "
            f"(expected ≥ {MIN_REPORT_BYTES:,}). Plotly.js may not be embedded."
        )

    def test_plotly_charts_present(self, cohort_output: Path) -> None:
        content = cohort_output.read_text(encoding="utf-8", errors="replace")
        assert "Plotly.newPlot" in content

    def test_no_external_cdn(self, cohort_output: Path) -> None:
        content = cohort_output.read_text(encoding="utf-8", errors="replace")
        assert "cdn.plot.ly" not in content, (
            "Report loads Plotly from CDN — it is not self-contained"
        )


@pytest.mark.integration
class TestPatientReports:
    def test_patient_html_in_process(self, integrated_result) -> None:
        patient_id = integrated_result.rows[0].patient
        html = render_patient_html(patient_id, integrated_result)
        assert "<html" in html.lower()
        assert len(html) >= MIN_PATIENT_BYTES
        assert "Plotly.newPlot" in html


@pytest.mark.integration
class TestErrorHandling:
    def test_missing_file_exits_nonzero(self, tmp_path: Path) -> None:
        result = subprocess.run(
            [
                sys.executable, "-m", "report_generator.generate_report",
                "--feature-table", str(DATA / "nonexistent.tsv"),
                "--taxonomy",      str(TAXONOMY),
                "--metadata",      str(METADATA),
                "--output",        str(tmp_path / "report.html"),
            ],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode != 0
