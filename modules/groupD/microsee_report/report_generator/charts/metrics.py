"""charts/metrics.py — alpha diversity metric labels and row-level value extraction.

Extracted from alpha.py so that stats.py can import metric helpers without
creating a circular dependency (stats ↔ alpha).  Both alpha.py and stats.py
now import from this module.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np

METRIC_LABELS: dict[str, str] = {
    "shannon":  "Shannon H′",
    "simpson":  "Simpson 1−D",
    "pielou":   "Pielou J′",
    "observed": "Observed Taxa",
    "faith_pd": "Faith PD",
}


def pielou_evenness(row: dict[str, Any], taxa: list[str]) -> float:
    """Compute Pielou J′ from family-level relative abundances stored in a row dict."""
    probs = np.array([float(row.get(t) or 0) / 100.0 for t in taxa])
    probs = probs[probs > 0]
    obs   = len(probs)
    if obs < 2:
        return 0.0
    return round(float(-np.sum(probs * np.log(probs + 1e-12))) / math.log(obs), 3)


def metric_value(row: dict[str, Any], metric: str, taxa: list[str] | None = None) -> float:
    """Read an alpha diversity metric value from a row dict.

    Falls back to computing pielou/observed from taxon abundances when the stored
    value is zero and taxa is provided.  This handles datasets where --alpha was
    not supplied and fallback estimates were used during integration.
    """
    v = float(row.get(metric) or 0)
    if v > 0:
        return v
    if taxa:
        if metric == "pielou":
            return pielou_evenness(row, taxa)
        if metric == "observed":
            return float(sum(1 for t in taxa if float(row.get(t) or 0) > 0.5))
    return 0.0
