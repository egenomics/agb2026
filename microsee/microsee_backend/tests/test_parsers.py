"""
tests/test_parsers.py

Unit tests for app/services/parsers.py.

Run with:  pytest tests/ -v
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import pytest
from app.services.parsers import (
    parse_feature_table,
    parse_taxonomy,
    parse_metadata,
    parse_alpha_diversity,
    parse_distance_matrix,
    integrate,
    _extract_family,
    _parse_time_days,
)


# ══════════════════════════════════════════════════════════════════════════════
# FIXTURES — minimal but realistic QIIME2-format strings
# ══════════════════════════════════════════════════════════════════════════════

FEATURE_TABLE_TSV = """\
# Constructed from biom file
#OTU ID\tEAA01_T0\tEAA01_T84\tWHY01_T0\tWHY01_T84
asv001\t3600\t3400\t4200\t4300
asv002\t2200\t2300\t2000\t1900
asv003\t1400\t1400\t1300\t1300
asv004\t900\t1000\t800\t800
asv005\t600\t700\t500\t500
asv006\t300\t300\t400\t400
"""

TAXONOMY_TSV = """\
Feature ID\tTaxon\tConfidence
asv001\td__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides\t0.99
asv002\td__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Lachnospira\t0.97
asv003\td__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Ruminococcus\t0.95
asv004\td__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella\t0.93
asv005\td__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium\t0.91
asv006\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia\t0.88
"""

METADATA_TSV = """\
sample-id\tgroup\ttimepoint\tpatient_id\tsixmwt\til18
EAA01_T0\tEAA\tT0\tEAA01\t310\t248
EAA01_T84\tEAA\tT84\tEAA01\t372\t192
WHY01_T0\tWhey\tT0\tWHY01\t315\t242
WHY01_T84\tWhey\tT84\tWHY01\t330\t240
"""

ALPHA_TSV = """\
sample-id\tshannon_entropy\tobserved_features\tfaith_pd\tpielou_evenness
EAA01_T0\t1.791\t6\t8.4\t0.998
EAA01_T84\t1.814\t6\t8.6\t0.999
WHY01_T0\t1.723\t5\t7.9\t0.982
WHY01_T84\t1.734\t5\t8.0\t0.983
"""

DISTANCE_MATRIX_TSV = """\
\tEAA01_T0\tEAA01_T84\tWHY01_T0\tWHY01_T84
EAA01_T0\t0.0\t0.05\t0.42\t0.41
EAA01_T84\t0.05\t0.0\t0.40\t0.39
WHY01_T0\t0.42\t0.40\t0.0\t0.04
WHY01_T84\t0.41\t0.39\t0.04\t0.0
"""


# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════

class TestHelpers:
    def test_extract_family_silva(self):
        t = "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides"
        assert _extract_family(t) == "Bacteroidaceae"

    def test_extract_family_greengenes(self):
        t = "k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae"
        assert _extract_family(t) == "Bacteroidaceae"

    def test_extract_family_unclassified(self):
        assert _extract_family("d__Bacteria;p__Firmicutes") == "Unclassified"

    def test_extract_family_unknown_value(self):
        assert _extract_family("d__Bacteria;f__uncultured") == "Unclassified"

    def test_extract_family_non_string(self):
        assert _extract_family(None) == "Unclassified"   # type: ignore

    def test_parse_time_t0(self):
        assert _parse_time_days("T0") == 0

    def test_parse_time_t84(self):
        assert _parse_time_days("T84") == 84

    def test_parse_time_week12(self):
        assert _parse_time_days("Week12") == 84

    def test_parse_time_numeric(self):
        assert _parse_time_days("0") == 0
        assert _parse_time_days("14") == 14

    def test_parse_time_baseline(self):
        assert _parse_time_days("baseline") == 0

    def test_parse_time_unparseable(self):
        assert _parse_time_days("visit_A") is None


# ══════════════════════════════════════════════════════════════════════════════
# FEATURE TABLE
# ══════════════════════════════════════════════════════════════════════════════

class TestParseFeatureTable:
    def test_basic_shape(self):
        result = parse_feature_table(FEATURE_TABLE_TSV)
        assert result.n_samples  == 4
        assert result.n_features == 6

    def test_sample_ids(self):
        result = parse_feature_table(FEATURE_TABLE_TSV)
        assert set(result.samples) == {"EAA01_T0", "EAA01_T84", "WHY01_T0", "WHY01_T84"}

    def test_feature_ids(self):
        result = parse_feature_table(FEATURE_TABLE_TSV)
        assert "asv001" in result.features

    def test_count_value(self):
        result = parse_feature_table(FEATURE_TABLE_TSV)
        assert result.counts["asv001"]["EAA01_T0"] == 3600.0

    def test_skips_biom_comment(self):
        # The # Constructed line should be silently ignored
        result = parse_feature_table(FEATURE_TABLE_TSV)
        assert result.n_features == 6  # not 7

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            parse_feature_table("")

    def test_malformed_raises(self):
        # pandas is very lenient — test a genuinely unparseable input
        with pytest.raises(ValueError):
            parse_feature_table("# Constructed from biom file\n")  # only comment, no data


# ══════════════════════════════════════════════════════════════════════════════
# TAXONOMY
# ══════════════════════════════════════════════════════════════════════════════

class TestParseTaxonomy:
    def test_family_extraction(self):
        result = parse_taxonomy(TAXONOMY_TSV)
        assert result.assignments["asv001"] == "Bacteroidaceae"
        assert result.assignments["asv002"] == "Lachnospiraceae"
        assert result.assignments["asv006"] == "Enterobacteriaceae"

    def test_unique_families(self):
        result = parse_taxonomy(TAXONOMY_TSV)
        assert "Bacteroidaceae" in result.unique_families
        assert "Ruminococcaceae" in result.unique_families

    def test_zero_unclassified(self):
        result = parse_taxonomy(TAXONOMY_TSV)
        assert result.unclassified_pct == 0.0

    def test_partial_unclassified(self):
        tsv = "Feature ID\tTaxon\nf1\td__Bacteria;p__Firmicutes\nf2\td__Bacteria;f__Bacteroidaceae"
        result = parse_taxonomy(tsv)
        assert result.unclassified_pct == 50.0

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            parse_taxonomy("")

    def test_missing_feature_id_column_raises(self):
        with pytest.raises(ValueError, match="Feature ID"):
            parse_taxonomy("Taxon\tConfidence\nd__Bacteria\t0.99")


# ══════════════════════════════════════════════════════════════════════════════
# METADATA
# ══════════════════════════════════════════════════════════════════════════════

class TestParseMetadata:
    def test_sample_count(self):
        result = parse_metadata(METADATA_TSV)
        assert result.n_samples == 4

    def test_group_derivation(self):
        result = parse_metadata(METADATA_TSV)
        by_id = {s.sample_id: s for s in result.samples}
        assert by_id["EAA01_T0"].group      == "EAA_T0"
        assert by_id["EAA01_T0"].base_group == "EAA"
        assert by_id["EAA01_T0"].timepoint  == "T0"

    def test_time_days(self):
        result = parse_metadata(METADATA_TSV)
        by_id  = {s.sample_id: s for s in result.samples}
        assert by_id["EAA01_T0"].time  == 0
        assert by_id["EAA01_T84"].time == 84

    def test_patient_id(self):
        result = parse_metadata(METADATA_TSV)
        by_id  = {s.sample_id: s for s in result.samples}
        assert by_id["EAA01_T0"].patient == "EAA01"

    def test_clinical_columns(self):
        result = parse_metadata(METADATA_TSV)
        assert result.has_clinical is True
        by_id  = {s.sample_id: s for s in result.samples}
        assert by_id["EAA01_T0"].sixmwt == 310.0
        assert by_id["EAA01_T0"].il18   == 248.0

    def test_groups_list(self):
        result = parse_metadata(METADATA_TSV)
        assert "EAA_T0"  in result.groups
        assert "Whey_T0" in result.groups

    def test_base_groups_list(self):
        result = parse_metadata(METADATA_TSV)
        assert "EAA"  in result.base_groups
        assert "Whey" in result.base_groups

    def test_missing_sample_id_raises(self):
        bad = "subject\tgroup\ns1\tEAA"
        with pytest.raises(ValueError, match="sample-id"):
            parse_metadata(bad)

    def test_skips_q2types_line(self):
        tsv_with_directive = "#q2:types\tcategorical\tcategorical\n" + METADATA_TSV
        result = parse_metadata(tsv_with_directive)
        assert result.n_samples == 4

    def test_no_clinical_columns(self):
        minimal = "sample-id\tgroup\ntimepoint\nEAA01_T0\tEAA\nT0"
        # Should parse without errors even without clinical columns
        tsv = "sample-id\tgroup\nEAA01_T0\tEAA\nWHY01_T0\tWhey"
        result = parse_metadata(tsv)
        assert result.has_clinical is False


# ══════════════════════════════════════════════════════════════════════════════
# ALPHA DIVERSITY
# ══════════════════════════════════════════════════════════════════════════════

class TestParseAlphaDiversity:
    def test_sample_count(self):
        result = parse_alpha_diversity(ALPHA_TSV)
        assert len(result.samples) == 4

    def test_shannon_values(self):
        result = parse_alpha_diversity(ALPHA_TSV)
        by_id  = {e.sample_id: e for e in result.samples}
        assert by_id["EAA01_T0"].shannon  == pytest.approx(1.791, abs=0.001)

    def test_metrics_detected(self):
        result = parse_alpha_diversity(ALPHA_TSV)
        assert "shannon"  in result.metrics_present
        assert "observed" in result.metrics_present
        assert "pielou"   in result.metrics_present
        assert "faith_pd" in result.metrics_present

    def test_single_metric_file(self):
        tsv = "sample-id\tshannon_entropy\nEAA01_T0\t1.791\nWHY01_T0\t1.723"
        result = parse_alpha_diversity(tsv)
        assert result.metrics_present == ["shannon"]
        assert len(result.samples) == 2

    def test_no_metric_columns_raises(self):
        tsv = "sample-id\tsome_other_column\nS1\t1.0"
        with pytest.raises(ValueError, match="No recognised"):
            parse_alpha_diversity(tsv)


# ══════════════════════════════════════════════════════════════════════════════
# DISTANCE MATRIX
# ══════════════════════════════════════════════════════════════════════════════

class TestParseDistanceMatrix:
    def test_shape(self):
        result = parse_distance_matrix(DISTANCE_MATRIX_TSV)
        assert result.n == 4
        assert len(result.samples) == 4
        assert len(result.matrix) == 4

    def test_sample_ids(self):
        result = parse_distance_matrix(DISTANCE_MATRIX_TSV)
        assert result.samples[0] == "EAA01_T0"

    def test_diagonal_zero(self):
        result = parse_distance_matrix(DISTANCE_MATRIX_TSV)
        for i in range(result.n):
            assert result.matrix[i][i] == pytest.approx(0.0)

    def test_symmetry(self):
        result = parse_distance_matrix(DISTANCE_MATRIX_TSV)
        for i in range(result.n):
            for j in range(result.n):
                assert result.matrix[i][j] == pytest.approx(result.matrix[j][i], abs=1e-4)

    def test_distance_value(self):
        result = parse_distance_matrix(DISTANCE_MATRIX_TSV)
        # EAA01_T0 vs EAA01_T84 should be ~0.05 (same patient)
        assert result.matrix[0][1] == pytest.approx(0.05)

    def test_too_few_rows_raises(self):
        with pytest.raises(ValueError, match="fewer than 2"):
            parse_distance_matrix("")

    def test_non_square_raises(self):
        # 2 column headers but 3 data rows = not square
        bad = "\tsample1\tsample2\nsample1\t0\t0.5\nsample2\t0.5\t0\nsample3\t0.3\t0.4"
        with pytest.raises(ValueError, match="square|rows|columns"):
            parse_distance_matrix(bad)


# ══════════════════════════════════════════════════════════════════════════════
# INTEGRATION
# ══════════════════════════════════════════════════════════════════════════════

class TestIntegrate:
    @pytest.fixture
    def parsed(self):
        feat  = parse_feature_table(FEATURE_TABLE_TSV)
        tax   = parse_taxonomy(TAXONOMY_TSV)
        meta  = parse_metadata(METADATA_TSV)
        alpha = parse_alpha_diversity(ALPHA_TSV)
        return feat, tax, meta, alpha

    def test_row_count(self, parsed):
        feat, tax, meta, alpha = parsed
        result = integrate(feat, tax, meta, alpha)
        assert result.n_samples == 4

    def test_taxa_present(self, parsed):
        feat, tax, meta, alpha = parsed
        result = integrate(feat, tax, meta, alpha)
        assert "Bacteroidaceae"   in result.taxa
        assert "Lachnospiraceae"  in result.taxa
        assert "Enterobacteriaceae" in result.taxa

    def test_taxa_sorted_by_mean_abundance(self, parsed):
        feat, tax, meta, alpha = parsed
        result = integrate(feat, tax, meta, alpha)
        # Bacteroidaceae has the highest counts so should be first
        assert result.taxa[0] == "Bacteroidaceae"

    def test_relative_abundance_sums_to_100(self, parsed):
        feat, tax, meta, alpha = parsed
        result = integrate(feat, tax, meta, alpha)
        row = result.rows[0]
        total = sum(
            getattr(row, t, 0) or row.model_extra.get(t, 0)
            for t in result.taxa
        )
        assert total == pytest.approx(100.0, abs=0.5)

    def test_alpha_diversity_joined(self, parsed):
        feat, tax, meta, alpha = parsed
        result = integrate(feat, tax, meta, alpha)
        by_id  = {r.sample_id: r for r in result.rows}
        assert by_id["EAA01_T0"].shannon == pytest.approx(1.791, abs=0.001)

    def test_clinical_data_joined(self, parsed):
        feat, tax, meta, alpha = parsed
        result = integrate(feat, tax, meta, alpha)
        by_id  = {r.sample_id: r for r in result.rows}
        assert by_id["EAA01_T0"].sixmwt == 310.0
        assert by_id["EAA01_T0"].il18   == 248.0

    def test_groups_correct(self, parsed):
        feat, tax, meta, alpha = parsed
        result = integrate(feat, tax, meta, alpha)
        assert "EAA_T0"  in result.groups
        assert "Whey_T0" in result.groups

    def test_has_clinical_true(self, parsed):
        feat, tax, meta, alpha = parsed
        result = integrate(feat, tax, meta, alpha)
        assert result.has_clinical is True

    def test_mismatched_sample_ids_raises(self):
        """Feature table with sample IDs that don't exist in metadata."""
        feat  = parse_feature_table(FEATURE_TABLE_TSV)
        tax   = parse_taxonomy(TAXONOMY_TSV)
        # Metadata with completely different sample IDs
        bad_meta = parse_metadata(
            "sample-id\tgroup\ntimepoint\nXXX_T0\tEAA\nXXX_T84\tEAA"
        )
        with pytest.raises(ValueError, match="No samples survived"):
            integrate(feat, tax, bad_meta)

    def test_without_alpha(self, parsed):
        """Shannon/Simpson should be computed from abundances when alpha is absent."""
        feat, tax, meta, _ = parsed
        result = integrate(feat, tax, meta, alpha=None)
        by_id  = {r.sample_id: r for r in result.rows}
        # Shannon should be > 0 (computed from relative abundances)
        assert by_id["EAA01_T0"].shannon > 0
        assert by_id["EAA01_T0"].simpson > 0
