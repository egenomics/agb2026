"""End-to-end tests — raw TSVs through the full pipeline to chart data and HTML.

Covers:
  - Full parse → integrate → compute_chart_data chain
  - Edge cases: single timepoint, single group, no alpha file
  - Validation: duplicate IDs, malformed input, all-zero samples
  - Scientific integrity: abundances, sorted taxa, BC bounds, JSON safety
  - Reproducibility: identical inputs → identical outputs
  - HTML rendering (integration-marked, needs bundled plotly.js)

Datasets
--------
E2E_*  : 2 groups × 2 patients × 2 timepoints (8 samples)  — main test dataset
SINGLE_TP_META  : only T0 rows  — graceful-degradation test
SINGLE_GRP_META : one treatment group — comparative-chart test
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from report_generator.charts import compute_chart_data, render_html
from report_generator.parsers import (
    integrate,
    parse_alpha_diversity,
    parse_feature_table,
    parse_metadata,
    parse_taxonomy,
)

# ── Shared inline datasets ────────────────────────────────────────────────────

E2E_FEATURE_TABLE = """\
#OTU ID\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8
ASV1\t100\t200\t 50\t 80\t300\t310\t290\t305
ASV2\t300\t150\t400\t200\t 80\t 70\t 90\t 75
ASV3\t 50\t 50\t250\t100\t200\t210\t190\t205
ASV4\t 20\t 80\t 30\t 60\t150\t160\t140\t155
"""

E2E_TAXONOMY = """\
Feature ID\tTaxon\tConfidence
ASV1\td__Bacteria;p__Firmicutes;f__Lachnospiraceae\t0.99
ASV2\td__Bacteria;p__Bacteroidota;f__Bacteroidaceae\t0.98
ASV3\td__Bacteria;p__Firmicutes;f__Ruminococcaceae\t0.97
ASV4\td__Bacteria;p__Proteobacteria;f__Enterobacteriaceae\t0.95
"""

E2E_METADATA = """\
sample-id\tsubject\tgroup\ttimepoint
S1\tPat1\tEAA\tT0
S2\tPat1\tEAA\tT84
S3\tPat2\tEAA\tT0
S4\tPat2\tEAA\tT84
S5\tPat3\tControl\tT0
S6\tPat3\tControl\tT84
S7\tPat4\tControl\tT0
S8\tPat4\tControl\tT84
"""

E2E_ALPHA = """\
sample-id\tshannon_entropy\tsimpson\tfaith_pd
S1\t1.5\t0.75\t8.2
S2\t1.8\t0.82\t9.1
S3\t1.3\t0.70\t7.8
S4\t1.6\t0.78\t8.5
S5\t2.1\t0.88\t11.3
S6\t2.0\t0.86\t10.9
S7\t1.9\t0.84\t10.5
S8\t2.2\t0.90\t11.8
"""

# Single-timepoint: only T0 rows, no T84
SINGLE_TP_META = """\
sample-id\tsubject\tgroup\ttimepoint
S1\tPat1\tEAA\tT0
S3\tPat2\tEAA\tT0
S5\tPat3\tControl\tT0
S7\tPat4\tControl\tT0
"""

# Single-group: all samples in one treatment group
SINGLE_GRP_META = """\
sample-id\tsubject\tgroup\ttimepoint
S1\tPat1\tEAA\tT0
S2\tPat1\tEAA\tT84
S3\tPat2\tEAA\tT0
S4\tPat2\tEAA\tT84
S5\tPat3\tEAA\tT0
S6\tPat3\tEAA\tT84
S7\tPat4\tEAA\tT0
S8\tPat4\tEAA\tT84
"""

# ── Fixtures ──────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def e2e_result():
    """Fully integrated 2-group result (2 patients/group, T0+T84, with alpha)."""
    feat = parse_feature_table(E2E_FEATURE_TABLE)
    tax = parse_taxonomy(E2E_TAXONOMY)
    meta = parse_metadata(E2E_METADATA)
    alpha = parse_alpha_diversity(E2E_ALPHA)
    return integrate(feat, tax, meta, alpha)


@pytest.fixture(scope="module")
def e2e_chart_data(e2e_result):
    return compute_chart_data(e2e_result)


# ── PRIORITY 1: Full pipeline ─────────────────────────────────────────────────


class TestFullPipeline:
    """Parse → integrate → compute_chart_data — all without HTML rendering."""

    def test_result_dimensions(self, e2e_result):
        assert e2e_result.n_samples == 8
        assert e2e_result.n_taxa == 4
        assert len(e2e_result.rows) == 8

    def test_groups_derived_correctly(self, e2e_result):
        base_groups = sorted({r.base_group for r in e2e_result.rows})
        assert base_groups == ["Control", "EAA"]
        full_groups = sorted({r.group for r in e2e_result.rows})
        assert "EAA_T0" in full_groups
        assert "EAA_T84" in full_groups
        assert "Control_T0" in full_groups

    def test_taxa_sorted_by_mean_abundance_descending(self, e2e_result):
        rows_d = [r.model_dump() for r in e2e_result.rows]
        means = {
            t: sum(float(r.get(t) or 0) for r in rows_d) / len(rows_d)
            for t in e2e_result.taxa
        }
        for i in range(len(e2e_result.taxa) - 1):
            assert means[e2e_result.taxa[i]] >= means[e2e_result.taxa[i + 1]], (
                f"Taxa not sorted by abundance: {e2e_result.taxa[i]} ({means[e2e_result.taxa[i]]:.2f}) "
                f"< {e2e_result.taxa[i + 1]} ({means[e2e_result.taxa[i + 1]]:.2f})"
            )

    def test_chart_data_required_sections_present(self, e2e_chart_data):
        required = [
            "meta",
            "taxonomy_views",
            "alpha_metrics",
            "pcoa_bray",
            "paired_slope",
            "stability_bar",
            "diff_abundance",
            "longitudinal",
            "stats_table",
            "permanova",
            "insights",
        ]
        missing = [k for k in required if k not in e2e_chart_data]
        assert not missing, f"Missing chart_data sections: {missing}"

    def test_meta_fields(self, e2e_chart_data):
        m = e2e_chart_data["meta"]
        assert m["n_samples"] == 8
        assert m["n_taxa"] == 4
        assert len(m["base_groups"]) == 2
        assert "has_clinical" in m
        assert isinstance(m["warnings"], list)

    def test_chart_data_json_serializable(self, e2e_chart_data):
        """json.dumps(..., allow_nan=False) must not raise — render_html depends on this."""
        try:
            json.dumps(e2e_chart_data, allow_nan=False)
        except (ValueError, TypeError) as exc:
            pytest.fail(
                f"chart_data contains NaN/Infinity or non-serializable value: {exc}\n"
                "This will crash render_html at the __DATA_JSON__ substitution step."
            )

    def test_insights_non_empty(self, e2e_chart_data):
        ins = e2e_chart_data["insights"]
        assert isinstance(ins, dict)
        assert "taxonomy" in ins, "Taxonomy insight missing"
        assert len(ins["taxonomy"]) > 10, "Taxonomy insight suspiciously short"


# ── PRIORITY 3: Reproducibility ───────────────────────────────────────────────


class TestReproducibility:
    """Same input must produce bit-identical chart data on repeated calls."""

    def test_chart_data_deterministic(self, e2e_result):
        d1 = json.dumps(compute_chart_data(e2e_result), sort_keys=True, allow_nan=False)
        d2 = json.dumps(compute_chart_data(e2e_result), sort_keys=True, allow_nan=False)
        assert d1 == d2, (
            "compute_chart_data produced different output on identical inputs. "
            "Look for non-seeded random calls (np.random without rng) or "
            "non-deterministic set/dict iteration."
        )

    def test_permanova_deterministic(self, e2e_result):
        cd1 = compute_chart_data(e2e_result)
        cd2 = compute_chart_data(e2e_result)
        p1 = [row[3] for row in cd1["permanova"]["rows"]]
        p2 = [row[3] for row in cd2["permanova"]["rows"]]
        assert p1 == p2, "PERMANOVA p-values differ between runs — RNG seed not deterministic"

    def test_taxa_order_stable(self, e2e_result):
        r1 = integrate(
            parse_feature_table(E2E_FEATURE_TABLE),
            parse_taxonomy(E2E_TAXONOMY),
            parse_metadata(E2E_METADATA),
            parse_alpha_diversity(E2E_ALPHA),
        )
        assert r1.taxa == e2e_result.taxa, "Taxa order changed between integrate() calls"


# ── PRIORITY 2: Edge cases — graceful degradation ─────────────────────────────


class TestEdgeCases:
    """Charts must not crash on unusual-but-valid inputs."""

    def test_no_alpha_file_uses_fallback(self):
        """Without --alpha, Shannon/Simpson are estimated from relative abundances."""
        feat = parse_feature_table(E2E_FEATURE_TABLE)
        tax = parse_taxonomy(E2E_TAXONOMY)
        meta = parse_metadata(E2E_METADATA)
        result = integrate(feat, tax, meta, alpha=None)
        assert result.n_samples == 8
        for row in result.rows:
            assert row.shannon >= 0, f"{row.sample_id} has negative Shannon"
            assert row.simpson >= 0, f"{row.sample_id} has negative Simpson"

    def test_single_timepoint_no_crash(self):
        """Dataset with only T0 — longitudinal/delta/stability charts must not raise."""
        feat = parse_feature_table(E2E_FEATURE_TABLE)
        tax = parse_taxonomy(E2E_TAXONOMY)
        meta = parse_metadata(SINGLE_TP_META)
        result = integrate(feat, tax, meta, alpha=None)
        assert result.n_samples == 4
        # must not raise — all charts degrade gracefully when there's no T84
        chart_data = compute_chart_data(result)
        json.dumps(chart_data, allow_nan=False)

    def test_single_group_no_crash(self):
        """Dataset with one treatment group — comparative/PERMANOVA charts must not raise."""
        feat = parse_feature_table(E2E_FEATURE_TABLE)
        tax = parse_taxonomy(E2E_TAXONOMY)
        meta = parse_metadata(SINGLE_GRP_META)
        result = integrate(feat, tax, meta, alpha=None)
        assert result.n_samples == 8
        chart_data = compute_chart_data(result)
        json.dumps(chart_data, allow_nan=False)

    def test_mismatched_sample_ids_produce_warning(self):
        """Samples in the feature table with no metadata row produce a warning."""
        # EXTRA_META omits S7 and S8
        extra_meta = """\
sample-id\tsubject\tgroup\ttimepoint
S1\tPat1\tEAA\tT0
S2\tPat1\tEAA\tT84
S3\tPat2\tEAA\tT0
S4\tPat2\tEAA\tT84
S5\tPat3\tControl\tT0
S6\tPat3\tControl\tT84
"""
        feat = parse_feature_table(E2E_FEATURE_TABLE)
        tax = parse_taxonomy(E2E_TAXONOMY)
        meta = parse_metadata(extra_meta)
        result = integrate(feat, tax, meta, alpha=None)
        assert any("metadata" in w.lower() or "no metadata" in w.lower() for w in result.warnings), (
            f"Expected a warning about missing metadata, got: {result.warnings}"
        )

    def test_unclassified_features_skipped(self):
        """Unclassified features are excluded from abundances without crashing."""
        tax_with_unknown = """\
Feature ID\tTaxon\tConfidence
ASV1\td__Bacteria;p__Firmicutes;f__Lachnospiraceae\t0.99
ASV2\tUnassigned\t0.10
ASV3\td__Bacteria;p__Firmicutes;f__Ruminococcaceae\t0.97
ASV4\td__Bacteria;p__Proteobacteria;f__Enterobacteriaceae\t0.95
"""
        feat = parse_feature_table(E2E_FEATURE_TABLE)
        tax = parse_taxonomy(tax_with_unknown)
        meta = parse_metadata(E2E_METADATA)
        result = integrate(feat, tax, meta, alpha=None)
        # ASV2 (Unclassified) must be excluded from taxa list
        assert "Unclassified" not in result.taxa
        assert result.n_taxa == 3  # only 3 classifiable families

    def test_utf8_bom_feature_table(self):
        """BOM-prefixed files (common from Windows/Excel exports) must parse correctly."""
        bom_table = "﻿" + E2E_FEATURE_TABLE
        result = parse_feature_table(bom_table)
        assert result.n_samples == 8

    def test_utf8_bom_metadata(self):
        bom_meta = "﻿" + E2E_METADATA
        result = parse_metadata(bom_meta)
        assert result.n_samples == 8


# ── PRIORITY 2: Validation ────────────────────────────────────────────────────


class TestValidation:
    """Invalid inputs must raise ValueError with a human-readable message."""

    def test_duplicate_sample_id_feature_table(self):
        dup = "#OTU ID\tS1\tS1\nASV1\t100\t200\n"
        with pytest.raises(ValueError, match="[Dd]uplicate"):
            parse_feature_table(dup)

    def test_duplicate_sample_id_metadata(self):
        dup = "sample-id\tgroup\nS1\tEAA\nS1\tEAA\n"
        with pytest.raises(ValueError, match="[Dd]uplicate"):
            parse_metadata(dup)

    def test_empty_feature_table_raises(self):
        with pytest.raises(ValueError):
            parse_feature_table("")

    def test_empty_taxonomy_raises(self):
        with pytest.raises(ValueError, match="[Vv]alid|[Rr]ow|[Ee]mpty"):
            parse_taxonomy("Feature ID\tTaxon\n")  # header only, no data

    def test_missing_taxonomy_column_raises(self):
        with pytest.raises(ValueError):
            parse_taxonomy("col1\tcol2\nA\tB\n")

    def test_empty_metadata_raises(self):
        with pytest.raises(ValueError):
            parse_metadata("sample-id\tgroup\n")  # header only

    def test_no_surviving_samples_raises(self):
        """When no feature-table sample matches any metadata sample, integrate must raise."""
        feat = parse_feature_table(E2E_FEATURE_TABLE)
        tax = parse_taxonomy(E2E_TAXONOMY)
        # Metadata with entirely different sample IDs (X1, X2) — none overlap with S1-S8
        bad_meta = parse_metadata("sample-id\tgroup\ttimepoint\nX1\tEAA\tT0\nX2\tEAA\tT84\n")
        with pytest.raises(ValueError, match="[Ss]ample"):
            integrate(feat, tax, bad_meta, alpha=None)

    def test_no_classifiable_features_raises(self):
        """When every feature is Unclassified, integrate must raise."""
        all_unknown = """\
Feature ID\tTaxon\tConfidence
ASV1\tUnassigned\t0.10
ASV2\tUnassigned\t0.10
ASV3\tUnassigned\t0.10
ASV4\tUnassigned\t0.10
"""
        feat = parse_feature_table(E2E_FEATURE_TABLE)
        tax = parse_taxonomy(all_unknown)
        meta = parse_metadata(E2E_METADATA)
        with pytest.raises(ValueError, match="[Ff]amil"):
            integrate(feat, tax, meta, alpha=None)


# ── PRIORITY 8: Scientific integrity ─────────────────────────────────────────


class TestScientificIntegrity:
    """Results must satisfy known microbiome-data invariants."""

    def test_relative_abundances_sum_to_100(self, e2e_result):
        rows_d = [r.model_dump() for r in e2e_result.rows]
        for r in rows_d:
            total = sum(float(r.get(t) or 0) for t in e2e_result.taxa)
            assert total == pytest.approx(100.0, abs=0.5), (
                f"{r['sample_id']} relative abundances sum to {total:.2f}%  (expected ≈100)"
            )

    def test_alpha_diversity_non_negative(self, e2e_result):
        for row in e2e_result.rows:
            assert row.shannon >= 0, f"{row.sample_id}: Shannon = {row.shannon}"
            assert row.simpson >= 0, f"{row.sample_id}: Simpson = {row.simpson}"
            assert row.simpson <= 1.0, f"{row.sample_id}: Simpson > 1 ({row.simpson})"

    def test_stability_bar_bc_in_unit_interval(self, e2e_chart_data):
        stab = e2e_chart_data.get("stability_bar", [])
        assert stab, "stability_bar is empty — expected data from 4 patients"
        for bc_val in stab[0].get("x", []):
            assert 0.0 <= float(bc_val) <= 1.0, (
                f"Bray-Curtis value out of [0, 1]: {bc_val}"
            )

    def test_permanova_r2_in_unit_interval(self, e2e_chart_data):
        pm = e2e_chart_data.get("permanova", {})
        for row in pm.get("rows", []):
            r2 = float(row[1])
            assert 0.0 <= r2 <= 1.0, f"PERMANOVA R² out of [0, 1]: {r2} (variable: {row[0]})"

    def test_pcoa_variance_explained_non_negative(self, e2e_chart_data):
        for key in ("pcoa_bray", "pcoa_jaccard"):
            d = e2e_chart_data.get(key, {})
            assert d.get("pct1", -1) >= 0.0, f"{key} pct1 is negative"
            assert d.get("pct2", -1) >= 0.0, f"{key} pct2 is negative"

    def test_stats_table_has_all_metrics(self, e2e_chart_data):
        from report_generator.charts.metrics import METRIC_LABELS
        st = e2e_chart_data["stats_table"]
        header_str = " ".join(str(h) for h in st["header"])
        for label in METRIC_LABELS.values():
            assert label.split()[0] in header_str, (
                f"Metric '{label}' missing from stats table header: {st['header']}"
            )

    def test_alpha_metrics_all_present(self, e2e_chart_data):
        am = e2e_chart_data["alpha_metrics"]
        for metric in ("shannon", "simpson", "pielou", "observed", "faith_pd"):
            assert metric in am, f"Alpha metric '{metric}' missing from chart_data"
            assert len(am[metric]["strip"]) > 0, f"No strip traces for metric '{metric}'"


# ── PRIORITY 1: HTML rendering (integration — needs plotly.js) ────────────────


@pytest.mark.integration
class TestHTMLOutput:
    """HTML output must be structurally complete and contain no unreplaced placeholders."""

    def test_cohort_html_no_placeholders(self, e2e_chart_data):
        html = render_html(e2e_chart_data)
        known_placeholders = [
            "__DATA_JSON__",
            "__INSIGHTS_JSON__",
            "__PATIENT_NAV__",
            "__CLINICAL_NAV__",
            "__WARNINGS__",
            "__N_SAMPLES__",
            "__GROUPS_STR__",
            "__FONT__",
            "__BG__",
        ]
        unreplaced = [p for p in known_placeholders if p in html]
        assert not unreplaced, f"Unreplaced template placeholders: {unreplaced}"

    def test_cohort_html_has_required_sections(self, e2e_chart_data):
        html = render_html(e2e_chart_data)
        for section_id in [
            "sec-taxonomy",
            "sec-alpha",
            "sec-beta",
            "sec-individual",
            "sec-comparative",
            "sec-longitudinal",
            "sec-stats",
        ]:
            assert section_id in html, f"Section '{section_id}' missing from HTML"

    def test_patient_nav_injected_when_mode_all(self, e2e_result):
        cd = compute_chart_data(e2e_result)
        cd["meta"]["mode"] = "all"
        cd["meta"]["patient_nav"] = [
            {"label": "Pat1", "href": "report_Pat1.html"},
            {"label": "Pat2", "href": "report_Pat2.html"},
        ]
        html = render_html(cd)
        # CSS always has `.nav-group-label{...}`; check the injected HTML element instead
        assert '<div class="nav-group-label">' in html, "Patient nav element not injected"
        assert "Pat1" in html
        assert "Pat2" in html

    def test_patient_nav_absent_when_mode_cohort(self, e2e_chart_data):
        html = render_html(e2e_chart_data)
        assert '<div class="nav-group-label">' not in html, (
            "Patient nav element should not appear in cohort-mode HTML"
        )

    def test_cohort_html_data_json_is_valid(self, e2e_chart_data):
        """The embedded __DATA_JSON__ must be valid JSON the browser can parse."""
        html = render_html(e2e_chart_data)
        # Extract the inline data JSON: template injects `var D = <json>;` on its own line.
        # Use newline as end delimiter — semicolons can appear inside JSON string values
        # (file paths, platform strings) and would break semicolon-based extraction.
        marker = "var D = "
        start = html.find(marker)
        assert start != -1, f"Could not find '{marker}' in rendered HTML"
        end = html.find("\n", start + len(marker))
        data_str = html[start + len(marker) : end if end != -1 else None].rstrip(";").strip()
        try:
            json.loads(data_str)
        except json.JSONDecodeError as exc:
            pytest.fail(f"Embedded data JSON is not valid: {exc}")


# ── Fixture-file smoke test ───────────────────────────────────────────────────

_DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def fixture_result():
    return integrate(
        parse_feature_table((_DATA / "feature-table.tsv").read_text(encoding="utf-8")),
        parse_taxonomy((_DATA / "taxonomy.tsv").read_text(encoding="utf-8")),
        parse_metadata((_DATA / "metadata.tsv").read_text(encoding="utf-8")),
        parse_alpha_diversity((_DATA / "alpha-diversity.tsv").read_text(encoding="utf-8")),
    )


class TestFixtureE2E:
    """Spot-checks on the realistic 12-patient fixture dataset."""

    def test_fixture_pipeline_dimensions(self, fixture_result):
        assert fixture_result.n_samples == 24
        assert fixture_result.n_taxa == 9
        assert fixture_result.has_clinical is True

    def test_fixture_chart_data_json_safe(self, fixture_result):
        cd = compute_chart_data(fixture_result)
        try:
            json.dumps(cd, allow_nan=False)
        except (ValueError, TypeError) as exc:
            pytest.fail(f"Fixture chart_data has NaN/Inf: {exc}")

    def test_fixture_deterministic(self, fixture_result):
        cd1 = json.dumps(compute_chart_data(fixture_result), sort_keys=True, allow_nan=False)
        cd2 = json.dumps(compute_chart_data(fixture_result), sort_keys=True, allow_nan=False)
        assert cd1 == cd2
