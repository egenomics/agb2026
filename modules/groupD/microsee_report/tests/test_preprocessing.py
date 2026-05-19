"""Unit tests for charts/preprocessing.py — shared data-preparation helpers."""


from report_generator.charts.preprocessing import (
    filter_rows,
    get_base_groups,
    get_patient_timepoints,
    get_unique_patients,
    sorted_timepoints,
)

# ---------------------------------------------------------------------------
# Minimal fixture rows — enough to exercise all helpers.
# Each dict mirrors the shape produced by parsers.integrate() → SampleRow.model_dump().
# ---------------------------------------------------------------------------

ROWS = [
    {"sample_id": "EAA01_T0",  "patient": "EAA01", "group": "EAA_T0",  "base_group": "EAA",  "timepoint": "T0",  "time": 0},
    {"sample_id": "EAA01_T84", "patient": "EAA01", "group": "EAA_T84", "base_group": "EAA",  "timepoint": "T84", "time": 84},
    {"sample_id": "EAA02_T0",  "patient": "EAA02", "group": "EAA_T0",  "base_group": "EAA",  "timepoint": "T0",  "time": 0},
    {"sample_id": "EAA02_T84", "patient": "EAA02", "group": "EAA_T84", "base_group": "EAA",  "timepoint": "T84", "time": 84},
    {"sample_id": "WHY01_T0",  "patient": "WHY01", "group": "Whey_T0", "base_group": "Whey", "timepoint": "T0",  "time": 0},
    {"sample_id": "WHY01_T84", "patient": "WHY01", "group": "Whey_T84","base_group": "Whey", "timepoint": "T84", "time": 84},
]


# ── sorted_timepoints ─────────────────────────────────────────────────────────

class TestSortedTimepoints:
    def test_basic_order(self):
        tps = sorted_timepoints(ROWS)
        assert tps[0] == "T0"
        assert tps[1] == "T84"

    def test_empty_rows(self):
        assert sorted_timepoints([]) == []

    def test_single_timepoint(self):
        rows = [{"timepoint": "T0", "time": 0}]
        assert sorted_timepoints(rows) == ["T0"]

    def test_missing_timepoint_skipped(self):
        rows = [
            {"timepoint": "T0", "time": 0},
            {"timepoint": None, "time": None},
            {"timepoint": "T42", "time": 42},
        ]
        tps = sorted_timepoints(rows)
        assert "T0" in tps
        assert "T42" in tps
        assert None not in tps

    def test_numeric_order_not_lexicographic(self):
        rows = [
            {"timepoint": "T10",  "time": 10},
            {"timepoint": "T2",   "time": 2},
            {"timepoint": "T100", "time": 100},
        ]
        tps = sorted_timepoints(rows)
        assert tps == ["T2", "T10", "T100"]


# ── get_base_groups ───────────────────────────────────────────────────────────

class TestGetBaseGroups:
    def test_returns_sorted_unique(self):
        bgs = get_base_groups(ROWS)
        assert bgs == ["EAA", "Whey"]

    def test_empty(self):
        assert get_base_groups([]) == []

    def test_single_group(self):
        rows = [{"base_group": "EAA"}, {"base_group": "EAA"}]
        assert get_base_groups(rows) == ["EAA"]

    def test_falls_back_to_group_key(self):
        rows = [{"group": "X"}, {"group": "Y"}]
        bgs = get_base_groups(rows)
        assert set(bgs) == {"X", "Y"}


# ── get_unique_patients ───────────────────────────────────────────────────────

class TestGetUniquePatients:
    def test_sorted_unique(self):
        patients = get_unique_patients(ROWS)
        assert patients == ["EAA01", "EAA02", "WHY01"]

    def test_empty(self):
        assert get_unique_patients([]) == []

    def test_single_patient_multiple_rows(self):
        rows = [{"patient": "P1"}, {"patient": "P1"}, {"patient": "P1"}]
        assert get_unique_patients(rows) == ["P1"]


# ── get_patient_timepoints ────────────────────────────────────────────────────

class TestGetPatientTimepoints:
    def test_returns_t0_and_t84(self):
        r0, r84 = get_patient_timepoints(ROWS, "EAA01")
        assert r0 is not None
        assert r84 is not None
        assert r0["timepoint"] == "T0"
        assert r84["timepoint"] == "T84"

    def test_missing_patient_returns_none_none(self):
        r0, r84 = get_patient_timepoints(ROWS, "NONEXISTENT")
        assert r0 is None
        assert r84 is None

    def test_only_t0_row(self):
        rows = [{"patient": "P1", "timepoint": "T0", "time": 0}]
        r0, r84 = get_patient_timepoints(rows, "P1")
        assert r0 is not None
        assert r84 is None

    def test_only_post_rows(self):
        rows = [{"patient": "P1", "timepoint": "T84", "time": 84}]
        r0, r84 = get_patient_timepoints(rows, "P1")
        assert r0 is None
        assert r84 is not None

    def test_multiple_post_rows_takes_latest(self):
        rows = [
            {"patient": "P1", "timepoint": "T0",  "time": 0},
            {"patient": "P1", "timepoint": "T42", "time": 42},
            {"patient": "P1", "timepoint": "T84", "time": 84},
        ]
        r0, r_latest = get_patient_timepoints(rows, "P1")
        assert r0 is not None and r_latest is not None
        assert r0["timepoint"] == "T0"
        assert r_latest["timepoint"] == "T84"

    def test_filters_by_patient(self):
        r0_a, _ = get_patient_timepoints(ROWS, "EAA01")
        r0_b, _ = get_patient_timepoints(ROWS, "EAA02")
        assert r0_a is not None and r0_b is not None
        assert r0_a["sample_id"] == "EAA01_T0"
        assert r0_b["sample_id"] == "EAA02_T0"


# ── filter_rows ───────────────────────────────────────────────────────────────

class TestFilterRows:
    def test_no_filters_returns_all(self):
        assert filter_rows(ROWS) == ROWS

    def test_filter_by_timepoint(self):
        filtered = filter_rows(ROWS, timepoint="T0")
        assert all(r["timepoint"] == "T0" for r in filtered)
        assert len(filtered) == 3  # EAA01, EAA02, WHY01

    def test_filter_by_base_group(self):
        filtered = filter_rows(ROWS, base_group="EAA")
        assert all(r["base_group"] == "EAA" for r in filtered)
        assert len(filtered) == 4

    def test_filter_both(self):
        filtered = filter_rows(ROWS, timepoint="T84", base_group="Whey")
        assert len(filtered) == 1
        assert filtered[0]["sample_id"] == "WHY01_T84"

    def test_no_match_returns_empty(self):
        assert filter_rows(ROWS, timepoint="T999") == []

    def test_empty_rows(self):
        assert filter_rows([], timepoint="T0") == []
