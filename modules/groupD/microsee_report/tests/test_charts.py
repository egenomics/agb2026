"""Smoke tests — verify chart builders return non-empty, well-shaped output.

Runs against both inline strings and realistic fixture TSV files.
"""
from pathlib import Path

import numpy as np
import pytest

from report_generator.parsers import (
    parse_feature_table, parse_taxonomy, parse_metadata,
    parse_alpha_diversity, integrate,
)
from report_generator.charts import compute_chart_data, render_html


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


@pytest.fixture(scope="module")
def chart_data():
    feat   = parse_feature_table(FEATURE_TABLE)
    tax    = parse_taxonomy(TAXONOMY)
    meta   = parse_metadata(METADATA)
    alpha  = parse_alpha_diversity(ALPHA)
    result = integrate(feat, tax, meta, alpha)
    return compute_chart_data(result)


def test_meta_fields(chart_data):
    m = chart_data["meta"]
    assert m["n_samples"] == 4
    assert m["n_taxa"] > 0
    assert len(m["groups"]) > 0


def test_taxonomy_views_present(chart_data):
    assert "both:all:all" in chart_data["taxonomy_views"]
    assert len(chart_data["taxonomy_views"]["both:all:all"]) > 0


def test_alpha_metrics_all_present(chart_data):
    for metric in ("shannon", "simpson", "pielou", "observed", "faith_pd"):
        assert metric in chart_data["alpha_metrics"]
        assert len(chart_data["alpha_metrics"][metric]["strip"]) > 0


def test_alpha_brackets_have_shapes(chart_data):
    bk = chart_data["alpha_metrics"]["shannon"]["brackets"]
    assert "shapes" in bk
    assert "annots" in bk
    assert "y_max" in bk


def test_pcoa_traces_present(chart_data):
    assert len(chart_data["pcoa_bray"]["traces"]) > 0
    assert chart_data["pcoa_bray"]["pct1"] >= 0


def test_lme_trajectory_all_metrics(chart_data):
    for metric in ("shannon", "simpson"):
        assert metric in chart_data["lme_trajectory"]
        assert len(chart_data["lme_trajectory"][metric]) > 0


def test_stats_table_all_metrics(chart_data):
    from report_generator.charts.alpha import METRIC_LABELS
    expected_cols = 2 + len(METRIC_LABELS)  # Group + N + one column per metric
    st = chart_data["stats_table"]
    assert len(st["header"]) == expected_cols
    assert len(st["rows"]) > 0
    assert all(len(row) == expected_cols for row in st["rows"])


def test_render_html_produces_valid_output(chart_data):
    html = render_html(chart_data)
    assert html.startswith("<!DOCTYPE html>")
    assert "Plotly" in html or "<script>" in html
    assert "__DATA_JSON__" not in html   # placeholder was replaced
    assert len(html) > 100_000          # should be a substantial file


def test_insights_injected(chart_data):
    """Dynamic insights must be present and non-empty after rendering."""
    assert "insights" in chart_data
    ins = chart_data["insights"]
    assert isinstance(ins, dict)
    assert "taxonomy" in ins
    assert len(ins["taxonomy"]) > 20


# ── Fixture-file tests (clinical-complete 4-patient dataset) ─────────────────

_DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def fixture_chart_data():
    """Chart data built from the realistic fixture TSV files (includes clinical)."""
    feat  = parse_feature_table((_DATA / "feature-table.tsv").read_text(encoding="utf-8"))
    tax   = parse_taxonomy((_DATA / "taxonomy.tsv").read_text(encoding="utf-8"))
    meta  = parse_metadata((_DATA / "metadata.tsv").read_text(encoding="utf-8"))
    alpha = parse_alpha_diversity((_DATA / "alpha-diversity.tsv").read_text(encoding="utf-8"))
    result = integrate(feat, tax, meta, alpha)
    return compute_chart_data(result)


class TestFixtureCharts:
    """Chart-builder smoke tests on the realistic 8-sample fixture dataset."""

    def test_meta_eight_samples(self, fixture_chart_data):
        assert fixture_chart_data["meta"]["n_samples"] == 8

    def test_clinical_sections_present(self, fixture_chart_data):
        assert "clinical_sixmwt" in fixture_chart_data
        assert "clinical_il18" in fixture_chart_data
        assert "taxa_clinical_heatmap" in fixture_chart_data

    def test_all_alpha_metrics_populated(self, fixture_chart_data):
        for metric in ("shannon", "simpson", "pielou", "observed", "faith_pd"):
            assert metric in fixture_chart_data["alpha_metrics"]
            strips = fixture_chart_data["alpha_metrics"][metric]["strip"]
            assert len(strips) > 0

    def test_stability_bar_all_patients(self, fixture_chart_data):
        stab = fixture_chart_data["stability_bar"]
        assert len(stab) > 0
        # 4 patients should each have a stability score
        assert len(stab[0]["y"]) == 4

    def test_insights_taxonomy_and_permanova(self, fixture_chart_data):
        ins = fixture_chart_data["insights"]
        assert "taxonomy" in ins and len(ins["taxonomy"]) > 10
        assert "permanova" in ins and len(ins["permanova"]) > 10

    def test_render_html_no_placeholders(self, fixture_chart_data):
        html = render_html(fixture_chart_data)
        assert "__DATA_JSON__" not in html
        assert "__INSIGHTS_JSON__" not in html
        assert len(html) > 100_000


# ── Distance matrix correctness ───────────────────────────────────────────────

from report_generator.charts.distances import bray_curtis_matrix, jaccard_matrix


class TestBrayCurtisMatrix:
    """Verify the vectorized bray_curtis_matrix against hand-computed values."""

    def test_known_value(self):
        # a=[3,1,2]  b=[2,4,0]
        # num = |3-2|+|1-4|+|2-0| = 1+3+2 = 6
        # den = (3+1+2)+(2+4+0)   = 6+6   = 12  →  BC = 0.5
        ab  = np.array([[3., 1., 2.], [2., 4., 0.]])
        mat = bray_curtis_matrix(ab)
        assert mat[0, 1] == pytest.approx(0.5)
        assert mat[1, 0] == pytest.approx(0.5)  # symmetric

    def test_identical_samples_zero(self):
        ab  = np.array([[1., 2., 3.], [1., 2., 3.]])
        mat = bray_curtis_matrix(ab)
        assert mat[0, 1] == pytest.approx(0.0)

    def test_completely_disjoint_one(self):
        # a=[1,0]  b=[0,1]: num=2, den=2  →  BC = 1.0
        ab  = np.array([[1., 0.], [0., 1.]])
        mat = bray_curtis_matrix(ab)
        assert mat[0, 1] == pytest.approx(1.0)

    def test_symmetric_diagonal_zero_bounded(self):
        rng = np.random.default_rng(42)
        ab  = rng.random((6, 10))
        mat = bray_curtis_matrix(ab)
        assert np.allclose(mat, mat.T), "not symmetric"
        assert np.allclose(np.diag(mat), 0.0), "diagonal not zero"
        assert np.all(mat >= 0) and np.all(mat <= 1), "values outside [0, 1]"


class TestJaccardMatrix:
    """Verify the vectorized jaccard_matrix against hand-computed values."""

    def test_completely_disjoint_one(self):
        # A=[100,0]  B=[0,100]: no shared presence  →  J = 1.0
        ab  = np.array([[100., 0.], [0., 100.]])
        mat = jaccard_matrix(ab, threshold=5.0)
        assert mat[0, 1] == pytest.approx(1.0)

    def test_identical_samples_zero(self):
        ab  = np.array([[80., 20., 0.], [80., 20., 0.]])
        mat = jaccard_matrix(ab, threshold=5.0)
        assert mat[0, 1] == pytest.approx(0.0)

    def test_partial_overlap(self):
        # A=[50,50, 0]  B=[50, 0,50]  (threshold=5)
        # pres_A=[T,T,F]  pres_B=[T,F,T]
        # both=1 (taxon 0)  union=3  →  J = 1 - 1/3 = 2/3
        ab  = np.array([[50., 50., 0.], [50., 0., 50.]])
        mat = jaccard_matrix(ab, threshold=5.0)
        assert mat[0, 1] == pytest.approx(2 / 3, abs=1e-6)

    def test_symmetric_diagonal_zero_bounded(self):
        rng = np.random.default_rng(99)
        ab  = rng.random((6, 10)) * 100
        mat = jaccard_matrix(ab, threshold=5.0)
        assert np.allclose(mat, mat.T), "not symmetric"
        assert np.allclose(np.diag(mat), 0.0), "diagonal not zero"
        assert np.all(mat >= 0) and np.all(mat <= 1), "values outside [0, 1]"
