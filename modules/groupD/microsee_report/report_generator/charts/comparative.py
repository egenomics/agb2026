"""charts/comparative.py — comparative chart builders (LFC, volcano, heatmap, correlation matrix)."""
from __future__ import annotations

import math
from typing import Any

import numpy as np

from .distances import rows_to_ab
from .stats_helpers import bh_fdr, welch_ttest_p, wilcoxon_p
from .utils import base_group_color, hex_rgba

_LFC_THRESH  = 0.5
_FDR_THRESH  = 0.1
_PSEUDO      = 0.1
_CLR_PSEUDO  = 0.01
_MIN_SAMPLES = 2


def _base_groups(rows: list[dict[str, Any]]) -> list[str]:
    return sorted({r.get("base_group", r["group"]) for r in rows})


def _clr(row: dict[str, Any], taxa: list[str]) -> dict[str, float]:
    """CLR-transform a single sample row (pseudocount 0.01)."""
    vals = {t: float(row.get(t) or 0) + _CLR_PSEUDO for t in taxa}
    gm   = math.exp(sum(math.log(v) for v in vals.values()) / len(vals))
    return {t: math.log(v / gm) for t, v in vals.items()}


def _ancom_color(diff: float, q: float, c: str) -> str:
    """Colour an ANCOM bar: significant up -> red, down -> blue, ns -> group colour."""
    if q < _FDR_THRESH and diff > 0:
        return "#D84E6A"
    if q < _FDR_THRESH and diff < 0:
        return "#4A7ED4"
    return hex_rgba(c, 0.45)


def build_diff_abundance(rows: list[dict[str, Any]], taxa: list[str]) -> list[dict[str, Any]]:
    """Log2 fold change per taxon: T84 vs T0 within each base group."""
    traces: list[dict[str, Any]] = []
    for bg in _base_groups(rows):
        t0_rows  = [
            r for r in rows
            if r.get("base_group", r["group"]) == bg and r.get("timepoint") == "T0"
        ]
        t84_rows = [
            r for r in rows
            if r.get("base_group", r["group"]) == bg and r.get("timepoint") == "T84"
        ]
        if not t0_rows or not t84_rows:
            continue
        stats = [
            (
                t,
                float(np.mean([float(r.get(t) or 0) for r in t0_rows])),
                float(np.mean([float(r.get(t) or 0) for r in t84_rows])),
            )
            for t in taxa
        ]
        lfc_vals = [
            round(math.log2((m84 + _PSEUDO) / (m0 + _PSEUDO)), 3)
            for _, m0, m84 in stats
        ]
        traces.append({
            "type": "bar", "orientation": "h", "name": bg,
            "x": lfc_vals,
            "y": [t for t, *_ in stats],
            "marker": {"color": ["#D84E6A" if v > 0 else "#4A7ED4" for v in lfc_vals]},
            "hovertemplate": "%{customdata}<extra></extra>",
            "customdata": [
                f"<b>{t}</b><br>T0: {m0:.2f}%  T84: {m84:.2f}%<br>Log2FC: {lfc_v:.3f}"
                for (t, m0, m84), lfc_v in zip(stats, lfc_vals, strict=False)
            ],
        })
    return traces


def build_volcano(rows: list[dict[str, Any]], taxa: list[str]) -> list[dict[str, Any]]:
    """Volcano: log2FC vs -log10(p); points coloured by BH-FDR q < 0.1."""
    bgs    = _base_groups(rows)
    traces: list[dict[str, Any]] = []
    for bg in bgs:
        t0_rows  = [
            r for r in rows
            if r.get("base_group", r["group"]) == bg and r.get("timepoint") == "T0"
        ]
        t84_rows = [
            r for r in rows
            if r.get("base_group", r["group"]) == bg and r.get("timepoint") == "T84"
        ]
        if not t0_rows or not t84_rows:
            continue
        c   = base_group_color(bg, bgs)
        raw: list[tuple[str, float, float]] = []
        for t in taxa:
            a_vals = [float(r.get(t) or 0) for r in t0_rows]
            b_vals = [float(r.get(t) or 0) for r in t84_rows]
            lfc    = round(
                math.log2(
                    (float(np.mean(b_vals)) + _PSEUDO)
                    / (float(np.mean(a_vals)) + _PSEUDO)
                ),
                3,
            )
            raw.append((t, lfc, welch_ttest_p(a_vals, b_vals)))
        p_vals = [p for _, _, p in raw]
        q_vals = bh_fdr(p_vals)
        vstats = [
            (t, lfc, round(-math.log10(max(p, 1e-10)), 3), q)
            for (t, lfc, p), q in zip(raw, q_vals, strict=False)
        ]
        sig = [abs(lfc) > _LFC_THRESH and q < _FDR_THRESH for _, lfc, _, q in vstats]
        traces.append({
            "type": "scatter", "mode": "markers+text", "name": bg,
            "x": [lfc for _, lfc, _, _ in vstats],
            "y": [nlp for _, _, nlp, _ in vstats],
            "text": [t for t, *_ in vstats],
            "textposition": "top center", "textfont": {"size": 9},
            "marker": {
                "color": ["#D84E6A" if s else hex_rgba(c, 0.5) for s in sig],
                "size": 10, "line": {"width": 1, "color": "white"},
            },
            "customdata": [[q] for _, _, _, q in vstats],
            "hovertemplate": (
                "<b>%{text}</b><br>Log2FC: %{x:.3f}<br>"
                "-log10(p): %{y:.2f}<br>"
                "FDR q: %{customdata[0]:.3f}<extra></extra>"
            ),
        })
    return traces


def build_ancom_style(rows: list[dict[str, Any]], taxa: list[str]) -> list[dict[str, Any]]:
    """CLR-transformed paired Wilcoxon (T0->T84) per taxon, per base group.

    Centered log-ratio is the compositionally-correct transform for differential
    abundance - equivalent in spirit to ANCOM/ANCOM-BC without the bias correction
    step.  Results displayed as an effect-size bar (CLR mean diff) coloured by FDR.
    """
    bgs    = _base_groups(rows)
    traces: list[dict[str, Any]] = []

    for bg in bgs:
        c        = base_group_color(bg, bgs)
        bg_rows  = [r for r in rows if r.get("base_group", r["group"]) == bg]
        clr_rows = {r["sample_id"]: _clr(r, taxa) for r in bg_rows}

        patients: list[str]                      = sorted({r["patient"] for r in bg_rows})
        taxon_stats: list[tuple[str, float, float]] = []
        p_vals_raw: list[float]                  = []

        for t in taxa:
            a: list[float] = []
            b: list[float] = []
            for pat in patients:
                t0r  = [
                    clr_rows[r["sample_id"]][t]
                    for r in bg_rows
                    if r["patient"] == pat and (r.get("time") or 0) == 0
                    and r["sample_id"] in clr_rows
                ]
                t84r = [
                    clr_rows[r["sample_id"]][t]
                    for r in bg_rows
                    if r["patient"] == pat and (r.get("time") or 0) > 0
                    and r["sample_id"] in clr_rows
                ]
                if t0r and t84r:
                    a.append(t0r[0])
                    b.append(t84r[0])
            diff  = round(float(np.mean(b)) - float(np.mean(a)), 4) if a else 0.0
            p_raw = wilcoxon_p(a, b)
            taxon_stats.append((t, diff, p_raw))
            p_vals_raw.append(p_raw)

        q_vals  = bh_fdr(p_vals_raw)
        results = [
            (t, diff, p, q)
            for (t, diff, p), q in zip(taxon_stats, q_vals, strict=False)
        ]
        traces.append({
            "type": "bar", "orientation": "h", "name": bg,
            "x": [diff for _, diff, _, _ in results],
            "y": [t for t, *_ in results],
            "marker": {"color": [_ancom_color(d, q, c) for _, d, _, q in results]},
            "customdata": [[p, q] for _, _, p, q in results],
            "hovertemplate": (
                "<b>%{y}</b><br>CLR mean diff: %{x:.3f}<br>"
                "p=%{customdata[0]:.3f}  FDR q=%{customdata[1]:.3f}<extra></extra>"
            ),
        })

    return traces


def build_heatmap(rows: list[dict[str, Any]], taxa: list[str]) -> list[dict[str, Any]]:
    """Raw abundance heatmap: samples x taxa, values as % relative abundance."""
    return [{
        "type": "heatmap",
        "x": list(taxa),
        "y": [r["sample_id"] for r in rows],
        "z": [[round(float(r.get(t) or 0), 2) for t in taxa] for r in rows],
        "colorscale": "Blues",
        "hovertemplate": "<b>%{y}</b> · <b>%{x}</b><br>%{z:.2f}%<extra></extra>",
    }]


def build_corr_matrix(rows: list[dict[str, Any]], taxa: list[str]) -> list[dict[str, Any]]:
    """Pearson correlation matrix across taxa."""
    ab = rows_to_ab(rows, taxa)
    if ab.shape[0] < _MIN_SAMPLES:
        return []
    corr = np.corrcoef(ab.T)
    return [{
        "type": "heatmap",
        "x": list(taxa), "y": list(taxa),
        "z": [
            [round(float(corr[i, j]), 3) for j in range(len(taxa))]
            for i in range(len(taxa))
        ],
        "colorscale": "RdBu", "zmid": 0, "zmin": -1, "zmax": 1,
        "hovertemplate": "<b>%{x}</b> vs <b>%{y}</b><br>r = %{z:.3f}<extra></extra>",
    }]
