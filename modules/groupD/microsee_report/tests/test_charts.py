"""Smoke tests — verify chart builders return non-empty, well-shaped output."""
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
S1\tPat1\tEAA_T0\tT0
S2\tPat1\tEAA_T84\tT84
S3\tPat2\tEAA_T0\tT0
S4\tPat2\tEAA_T84\tT84
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
    st = chart_data["stats_table"]
    assert len(st["header"]) == 7   # Group + N + 5 metrics
    assert len(st["rows"]) > 0
    assert all(len(row) == 7 for row in st["rows"])


def test_render_html_produces_valid_output(chart_data):
    html = render_html(chart_data)
    assert html.startswith("<!DOCTYPE html>")
    assert "Plotly" in html or "<script>" in html
    assert "__DATA_JSON__" not in html   # placeholder was replaced
    assert len(html) > 100_000          # should be a substantial file
