"""Unit tests for report_generator parsers."""
import pytest
from report_generator.parsers import (
    parse_feature_table,
    parse_taxonomy,
    parse_metadata,
    parse_alpha_diversity,
    integrate,
)


FEATURE_TABLE = """\
#OTU ID\tS1\tS2\tS3
ASV1\t100\t200\t50
ASV2\t300\t150\t400
ASV3\t0\t50\t250
"""

TAXONOMY = """\
Feature ID\tTaxon\tConfidence
ASV1\td__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lachnospiraceae\t0.99
ASV2\td__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae\t0.98
ASV3\td__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Ruminococcaceae\t0.97
"""

METADATA = """\
sample-id\tsubject\tgroup\ttimepoint
S1\tPat1\tEAA\tT0
S2\tPat1\tEAA\tT84
S3\tPat2\tEAA\tT0
"""

ALPHA = """\
sample-id\tshannon_entropy\tsimpson
S1\t1.5\t0.75
S2\t1.8\t0.82
S3\t1.3\t0.70
"""


def test_parse_feature_table_basic():
    result = parse_feature_table(FEATURE_TABLE)
    assert set(result.samples) == {"S1", "S2", "S3"}
    assert len(result.features) == 3
    assert result.n_samples == 3
    assert result.n_features == 3


def test_parse_feature_table_counts():
    result = parse_feature_table(FEATURE_TABLE)
    assert result.counts["ASV1"]["S1"] == 100
    assert result.counts["ASV2"]["S3"] == 400


def test_parse_taxonomy_family_level():
    result = parse_taxonomy(TAXONOMY)
    assert "Lachnospiraceae" in result.assignments.values()
    assert "Bacteroidaceae" in result.assignments.values()
    assert result.unclassified_pct == 0.0


def test_parse_metadata_groups():
    result = parse_metadata(METADATA)
    assert result.n_samples == 3
    assert any("EAA" in g for g in result.groups)
    assert any("T0" in g for g in result.groups)
    assert any("T84" in g for g in result.groups)
    assert result.has_clinical is False


def test_parse_alpha_diversity():
    result = parse_alpha_diversity(ALPHA)
    assert len(result.samples) == 3
    entry = next(e for e in result.samples if e.sample_id == "S1")
    assert entry.shannon == pytest.approx(1.5)
    assert entry.simpson == pytest.approx(0.75)


def test_integrate_produces_sample_rows():
    feat  = parse_feature_table(FEATURE_TABLE)
    tax   = parse_taxonomy(TAXONOMY)
    meta  = parse_metadata(METADATA)
    alpha = parse_alpha_diversity(ALPHA)
    result = integrate(feat, tax, meta, alpha)

    assert result.n_samples == 3
    assert result.n_taxa > 0
    assert len(result.rows) == 3
    assert result.has_clinical is False


def test_integrate_relative_abundances_sum_to_100():
    feat   = parse_feature_table(FEATURE_TABLE)
    tax    = parse_taxonomy(TAXONOMY)
    meta   = parse_metadata(METADATA)
    result = integrate(feat, tax, meta, None)

    for row in result.rows:
        total = sum(float(row.model_dump().get(t, 0)) for t in result.taxa)
        assert total == pytest.approx(100.0, abs=0.5), f"{row.sample_id} sums to {total}"


def test_empty_feature_table_raises():
    with pytest.raises(ValueError):
        parse_feature_table("")


def test_missing_taxonomy_column_raises():
    with pytest.raises(ValueError):
        parse_taxonomy("col1\tcol2\nA\tB\n")
