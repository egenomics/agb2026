"""charts/preprocessing.py — shared data-preparation helpers for chart builders.

All functions operate on the canonical list[dict[str, Any]] row representation produced by
parsers.integrate().  Centralising these avoids duplicating the same filtering
and pairing logic across beta, individual, comparative, clinical, stats, and
renderer — where the T0/T84 patient-pair pattern previously appeared eight times.
"""
from __future__ import annotations

from typing import Any

# ── Timepoint utilities ───────────────────────────────────────────────────────

def sorted_timepoints(rows: list[dict[str, Any]]) -> list[str]:
    """Return unique timepoints sorted by their numeric `time` field."""
    tp_time: dict[str, int] = {}
    for r in rows:
        tp = r.get("timepoint")
        if tp:
            tp_time[tp] = int(r.get("time") or 0)
    return sorted(tp_time.keys(), key=lambda tp: tp_time[tp])


# ── Group utilities ───────────────────────────────────────────────────────────

def get_base_groups(rows: list[dict[str, Any]]) -> list[str]:
    """Return sorted unique base_group values from a row list."""
    return sorted(set(r.get("base_group", r.get("group", "")) for r in rows))


def get_unique_patients(rows: list[dict[str, Any]]) -> list[str]:
    """Return sorted unique patient IDs."""
    return sorted(set(r["patient"] for r in rows))


# ── Patient-level time pairing ────────────────────────────────────────────────

def get_patient_timepoints(
    rows: list[dict[str, Any]], patient: str
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    """Return (T0 row, latest post-T0 row) for a patient.

    Returns None for a position when no matching row exists.  Callers should
    check for None before using the results, e.g.::

        r0, r84 = get_patient_timepoints(rows, patient)
        if r0 is None or r84 is None:
            continue
    """
    t0_list  = [r for r in rows if r["patient"] == patient and (r.get("time") or 0) == 0]
    t84_list = sorted(
        [r for r in rows if r["patient"] == patient and (r.get("time") or 0) > 0],
        key=lambda r: r.get("time") or 0,
        reverse=True,
    )
    return (t0_list[0] if t0_list else None), (t84_list[0] if t84_list else None)


# ── Row filtering ─────────────────────────────────────────────────────────────

def filter_rows(
    rows: list[dict[str, Any]],
    *,
    timepoint: str | None = None,
    base_group: str | None = None,
) -> list[dict[str, Any]]:
    """Filter a row list by timepoint and/or base_group.

    Both filters are ANDed when both are provided.  Returns the original list
    unmodified when neither filter is specified.
    """
    result = rows
    if timepoint:
        result = [r for r in result if r.get("timepoint") == timepoint]
    if base_group:
        result = [r for r in result if r.get("base_group", r.get("group")) == base_group]
    return result
