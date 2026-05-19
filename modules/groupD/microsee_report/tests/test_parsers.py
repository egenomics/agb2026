"""Unit tests for report_generator parsers — inline strings + fixture files."""
from pathlib import Path

import pytest

from report_generator.parsers import (
    integrate,
    parse_alpha_diversity,
    parse_feature_table,
    parse_metadata,
    parse_taxonomy,
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


# ── Fixture-file tests (realistic 12-patient x 2-timepoint dataset) ───────────

_DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def file_feat() -> str:
    """Read feature table from fixture file."""
    return (_DATA / "feature-table.tsv").read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def file_tax() -> str:
    """Read taxonomy from fixture file."""
    return (_DATA / "taxonomy.tsv").read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def file_meta() -> str:
    """Read metadata (with clinical columns) from fixture file."""
    return (_DATA / "metadata.tsv").read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def file_alpha() -> str:
    """Read alpha diversity (all 5 metrics) from fixture file."""
    return (_DATA / "alpha-diversity.tsv").read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def file_result(file_feat, file_tax, file_meta, file_alpha):
    """Fully integrated result built from fixture files."""
    return integrate(
        parse_feature_table(file_feat),
        parse_taxonomy(file_tax),
        parse_metadata(file_meta),
        parse_alpha_diversity(file_alpha),
    )


class TestFixtureFiles:
    """Replicate core parser assertions using realistic fixture TSV files."""

    def test_feature_table_dimensions(self, file_feat):
        result = parse_feature_table(file_feat)
        assert result.n_samples == 24
        assert result.n_features == 9

    def test_feature_table_sample_ids(self, file_feat):
        result = parse_feature_table(file_feat)
        assert "EAA01_T0" in result.samples
        assert "WHY06_T84" in result.samples

    def test_taxonomy_nine_families(self, file_tax):
        result = parse_taxonomy(file_tax)
        assert len(result.assignments) == 9
        assert "Lachnospiraceae" in result.assignments.values()
        assert "Akkermansiaceae" in result.assignments.values()
        assert "Oscillospiraceae" in result.assignments.values()
        assert result.unclassified_pct == 0.0

    def test_metadata_clinical_detected(self, file_meta):
        result = parse_metadata(file_meta)
        assert result.n_samples == 24
        assert result.has_clinical is True  # sixmwt + il18 columns present

    def test_metadata_groups_present(self, file_meta):
        result = parse_metadata(file_meta)
        names = " ".join(result.groups)
        assert "EAA" in names
        assert "Whey" in names

    def test_alpha_all_five_metrics(self, file_alpha):
        result = parse_alpha_diversity(file_alpha)
        assert len(result.samples) == 24
        s = next(e for e in result.samples if e.sample_id == "EAA01_T0")
        assert s.shannon == pytest.approx(1.791)
        assert s.simpson == pytest.approx(0.786)
        assert s.faith_pd == pytest.approx(9.2)

    def test_integrate_clinical_flag(self, file_result):
        assert file_result.has_clinical is True

    def test_integrate_twenty_four_rows(self, file_result):
        assert file_result.n_samples == 24
        assert len(file_result.rows) == 24

    def test_integrate_abundances_sum_to_100(self, file_result):
        for row in file_result.rows:
            total = sum(float(row.model_dump().get(t, 0)) for t in file_result.taxa)
            assert total == pytest.approx(100.0, abs=0.5), \
                f"{row.sample_id} abundances sum to {total:.2f}"
