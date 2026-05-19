"""charts/individual.py — per-patient chart builders."""
from __future__ import annotations

import math
from typing import Any

import numpy as np

from .distances import bray_curtis_matrix, pcoa, rows_to_ab
from .preprocessing import (
    get_base_groups,
    get_patient_timepoints,
    get_unique_patients,
    sorted_timepoints,
)
from .utils import base_group_color, hex_rgba, taxon_color


# ── Paired Slopegraph ─────────────────────────────────────────────────────────

def build_paired_slope(rows: list[dict[str, Any]], groups: list[str]) -> dict[str, list[dict[str, Any]]]:
    """Per-patient slopegraph for Shannon and Simpson — returns dict keyed by metric."""
    return {m: _slope_traces(rows, m) for m in ("shannon", "simpson")}


def _slope_traces(rows: list[dict[str, Any]], metric: str) -> list[dict[str, Any]]:
    timepoints  = sorted_timepoints(rows)
    base_groups = get_base_groups(rows)
    traces: list[dict[str, Any]] = []
    for p in get_unique_patients(rows):
        p_rows = sorted([r for r in rows if r["patient"] == p],
                        key=lambda r: r.get("time") or 0)
        if len(p_rows) < 2:
            continue
        bg = p_rows[0].get("base_group", p_rows[0]["group"])
        c  = base_group_color(bg, base_groups)
        traces.append({
            "type": "scatter", "mode": "lines+markers", "name": p,
            "x": [r["timepoint"] for r in p_rows],
            "y": [float(r.get(metric) or 0) for r in p_rows],
            "line": {"color": hex_rgba(c, 0.55), "width": 1.5},
            "marker": {"size": 7, "color": c},
            "showlegend": False,
            "hovertemplate": f"<b>{p}</b><br>%{{x}}<br>{metric}: %{{y:.3f}}<extra></extra>",
        })
    for bg in base_groups:
        c = base_group_color(bg, base_groups)
        means:   list[float] = []
        xs_used: list[str]   = []
        for tp in timepoints:
            vals = [float(r.get(metric) or 0) for r in rows
                    if r.get("base_group", r["group"]) == bg and r["timepoint"] == tp]
            if vals:
                means.append(round(float(np.mean(vals)), 4))
                xs_used.append(tp)
        if means:
            traces.append({
                "type": "scatter", "mode": "lines+markers",
                "name": f"{bg} mean", "x": xs_used, "y": means,
                "line": {"color": c, "width": 3, "dash": "dash"},
                "marker": {"size": 10, "color": c, "symbol": "diamond"},
                "hovertemplate": f"<b>{bg} mean</b><br>%{{x}}<br>{metric}: %{{y:.3f}}<extra></extra>",
            })
    return traces


# ── Stability Bar ─────────────────────────────────────────────────────────────

def build_stability_bar(rows: list[dict[str, Any]], taxa: list[str]) -> list[dict[str, Any]]:
    """Per-patient Bray-Curtis dissimilarity T0 vs T84, sorted horizontal bar."""
    base_groups = get_base_groups(rows)
    scores: list[dict[str, Any]] = []
    for p in get_unique_patients(rows):
        r0, r84 = get_patient_timepoints(rows, p)
        if r0 is None or r84 is None:
            continue
        ab0  = np.array([float(r0.get(t) or 0) for t in taxa])
        ab84 = np.array([float(r84.get(t) or 0) for t in taxa])
        num  = float(np.sum(np.abs(ab0 - ab84)))
        den  = float(np.sum(ab0 + ab84))
        bc   = round(num / den, 3) if den > 0 else 0.0
        bg   = r0.get("base_group", r0.get("group", ""))
        scores.append({"patient": p, "bc": bc, "base_group": bg})
    scores.sort(key=lambda s: s["bc"])
    return [{
        "type": "bar", "orientation": "h",
        "x": [s["bc"] for s in scores],
        "y": [s["patient"] for s in scores],
        "marker": {"color": [hex_rgba(base_group_color(s["base_group"], base_groups), 0.8) for s in scores]},
        "text": [str(s["bc"]) for s in scores],
        "textposition": "outside",
        "hovertemplate": "<b>%{y}</b><br>BC dissimilarity: %{x:.3f}<extra></extra>",
    }]


# ── Diversity Rank ────────────────────────────────────────────────────────────

def build_diversity_rank(rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    """Samples ranked low→high for Shannon and Simpson — returns dict keyed by metric."""
    return {m: _rank_traces(rows, m) for m in ("shannon", "simpson")}


def _rank_traces(rows: list[dict[str, Any]], metric: str) -> list[dict[str, Any]]:
    base_groups = get_base_groups(rows)
    sorted_rows = sorted(rows, key=lambda r: float(r.get(metric) or 0))
    traces: list[dict[str, Any]] = []
    for bg in base_groups:
        c   = base_group_color(bg, base_groups)
        pts = [(i, r) for i, r in enumerate(sorted_rows)
               if r.get("base_group", r.get("group")) == bg]
        t0  = [(i, r) for i, r in pts if (r.get("time") or 0) == 0]
        t84 = [(i, r) for i, r in pts if (r.get("time") or 0) > 0]
        if t0:
            traces.append({
                "type": "scatter", "mode": "markers", "name": f"{bg} T0",
                "x": [i for i, _ in t0],
                "y": [float(r.get(metric) or 0) for _, r in t0],
                "text": [r["sample_id"] for _, r in t0],
                "marker": {"color": hex_rgba(c, 0.6), "size": 8, "symbol": "circle",
                           "line": {"width": 1, "color": "white"}},
                "hovertemplate": "<b>%{text}</b>: %{y:.3f}<extra></extra>",
            })
        if t84:
            traces.append({
                "type": "scatter", "mode": "markers", "name": f"{bg} T84",
                "x": [i for i, _ in t84],
                "y": [float(r.get(metric) or 0) for _, r in t84],
                "text": [r["sample_id"] for _, r in t84],
                "marker": {"color": c, "size": 9, "symbol": "diamond",
                           "line": {"width": 1, "color": "white"}},
                "hovertemplate": "<b>%{text}</b>: %{y:.3f}<extra></extra>",
            })
    return traces


# ── Patient Radar Profiles ────────────────────────────────────────────────────

def build_patient_radar_profiles(rows: list[dict[str, Any]], taxa: list[str]) -> dict[str, Any]:
    """Pre-compute per-patient radar traces + composition table for all patients.

    Returns::
        {
          "patients": [str, ...],           # ordered patient list
          "profiles": {
            "<patient_id>": {
              "traces": [Plotly trace, ...], # T0 filled, T84 dashed, group mean dotted
              "table":  [{"family", "v0", "v84", "delta"}, ...],
              "group":  str,
            }
          }
        }
    """
    base_groups = get_base_groups(rows)
    short_taxa  = [t.replace("aceae", "")[:12] for t in taxa]

    # Pre-compute group mean T0 as the dotted reference line
    gmeans: dict[str, list[float]] = {}
    for bg in base_groups:
        t0s = [r for r in rows
               if r.get("base_group", r.get("group")) == bg and (r.get("time") or 0) == 0]
        if t0s:
            gmeans[bg] = [
                round(float(np.mean(
                    [float(r.get(t) or 0) / (sum(float(r.get(tt) or 0) for tt in taxa) or 1) * 100
                     for r in t0s]
                )), 1)
                for t in taxa
            ]

    patients = get_unique_patients(rows)
    profiles: dict[str, dict[str, Any]] = {}

    for p in patients:
        r0, r84 = get_patient_timepoints(rows, p)
        if r0 is None:
            continue
        bg = r0.get("base_group", r0.get("group", ""))
        c  = base_group_color(bg, base_groups)

        tot0 = sum(float(r0.get(t) or 0) for t in taxa) or 1.0
        v0   = [round(float(r0.get(t) or 0) / tot0 * 100, 1) for t in taxa]

        v84: list[float] | None = None
        if r84 is not None:
            tot84 = sum(float(r84.get(t) or 0) for t in taxa) or 1.0
            v84   = [round(float(r84.get(t) or 0) / tot84 * 100, 1) for t in taxa]

        traces: list[dict[str, Any]] = [{
            "type": "scatterpolar", "fill": "toself",
            "r": v0 + [v0[0]], "theta": short_taxa + [short_taxa[0]],
            "fillcolor": hex_rgba(c, 0.2),
            "line": {"color": hex_rgba(c, 0.8), "width": 2},
            "name": f"{p} T0",
            "hovertemplate": "<b>%{theta}</b>: %{r:.1f}%<extra></extra>",
        }]
        if v84 is not None:
            traces.append({
                "type": "scatterpolar", "fill": "none",
                "r": v84 + [v84[0]], "theta": short_taxa + [short_taxa[0]],
                "line": {"color": c, "width": 2.5, "dash": "dash"},
                "name": f"{p} T84",
                "hovertemplate": "<b>%{theta}</b>: %{r:.1f}%<extra></extra>",
            })
        if bg in gmeans:
            gm = gmeans[bg]
            traces.append({
                "type": "scatterpolar", "fill": "none",
                "r": gm + [gm[0]], "theta": short_taxa + [short_taxa[0]],
                "line": {"color": "#C4A0A8", "width": 1.5, "dash": "dot"},
                "name": f"{bg} mean T0",
                "hovertemplate": "<b>%{theta}</b>: %{r:.1f}%<extra></extra>",
            })

        table: list[dict[str, Any]] = []
        for i, t in enumerate(taxa):
            v84_i: float | None = v84[i] if v84 is not None else None
            delta: float | None = round(v84_i - v0[i], 1) if v84_i is not None else None
            table.append({
                "family": t.replace("aceae", ""),
                "v0": v0[i],
                "v84": v84_i,
                "delta": delta,
            })
        profiles[p] = {"traces": traces, "table": table, "group": bg}

    return {"patients": patients, "profiles": profiles}


# ── Faceted small multiples ───────────────────────────────────────────────────

def build_faceted_composition(rows: list[dict[str, Any]], taxa: list[str]) -> dict[str, Any]:
    """True small-multiples: one subplot per patient, T0 vs T84 stacked bars."""
    patients = get_unique_patients(rows)
    n        = len(patients)
    n_cols   = min(4, n)
    n_rows   = math.ceil(n / n_cols)
    col_gap, row_gap = 0.04, 0.10
    col_w = (1.0 - col_gap * (n_cols - 1)) / n_cols
    row_h = (1.0 - row_gap * (n_rows - 1)) / n_rows

    data: list[dict[str, Any]] = []
    annotations: list[dict[str, Any]] = []
    layout: dict[str, Any] = {
        "barmode": "stack", "showlegend": True,
        "height": max(360, n_rows * 240),
        "margin": {"l": 40, "r": 10, "t": 20, "b": 40},
    }
    shown: set[str] = set()

    for i, p in enumerate(patients):
        ri, ci  = divmod(i, n_cols)
        ax_suf  = "" if i == 0 else str(i + 1)
        x0 = round(ci * (col_w + col_gap), 4)
        x1 = round(x0 + col_w, 4)
        y0 = round((n_rows - 1 - ri) * (row_h + row_gap), 4)
        y1 = round(y0 + row_h, 4)

        layout[f"xaxis{ax_suf}"] = {
            "domain": [x0, x1], "anchor": f"y{ax_suf}",
            "tickfont": {"size": 9}, "tickangle": 0,
        }
        layout[f"yaxis{ax_suf}"] = {
            "domain": [y0, y1], "anchor": f"x{ax_suf}",
            "range": [0, 100], "ticksuffix": "%", "tickfont": {"size": 8},
        }
        annotations.append({
            "text": f"<b>{p}</b>",
            "x": (x0 + x1) / 2, "y": y1 + 0.005,
            "xref": "paper", "yref": "paper",
            "showarrow": False, "xanchor": "center", "yanchor": "bottom",
            "font": {"size": 10},
        })

        r0, r84 = get_patient_timepoints(rows, p)
        for r, tp in [(r0, "T0"), (r84, "T84")]:
            if not r:
                continue
            tot = sum(float(r.get(t) or 0) for t in taxa) or 1.0
            for t in taxa:
                pct = round(float(r.get(t) or 0) / tot * 100, 1)
                if pct == 0:
                    continue
                c = taxon_color(t)
                data.append({
                    "type": "bar", "name": t,
                    "x": [tp], "y": [pct],
                    "xaxis": f"x{ax_suf}", "yaxis": f"y{ax_suf}",
                    "marker": {"color": hex_rgba(c, 0.9 if tp == "T84" else 0.65)},
                    "legendgroup": t, "showlegend": t not in shown,
                    "hovertemplate": f"<b>{p} {tp}</b><br>{t}: %{{y:.1f}}%<extra></extra>",
                })
                shown.add(t)

    layout["annotations"] = annotations
    return {"data": data, "layout": layout}


# ── NMDS Trajectories ─────────────────────────────────────────────────────────

def build_nmds_trajectories(
    rows: list[dict[str, Any]],
    taxa: list[str],
    bc_matrix: np.ndarray | None = None,
) -> tuple[list[dict[str, Any]], float, float]:
    """PCoA (Bray-Curtis) with T0→T84 arrows per patient."""
    base_groups = get_base_groups(rows)
    mat = bc_matrix if bc_matrix is not None else bray_curtis_matrix(rows_to_ab(rows, taxa))
    xs, ys, pct1, pct2 = pcoa(mat)
    coord_map = {r["sample_id"]: (xs[i], ys[i]) for i, r in enumerate(rows)}
    traces: list[dict[str, Any]] = []
    for p in get_unique_patients(rows):
        r0, r84 = get_patient_timepoints(rows, p)
        if r0 is None or r84 is None:
            continue
        c0 = coord_map.get(r0["sample_id"])
        c1 = coord_map.get(r84["sample_id"])
        if not c0 or not c1:
            continue
        bg = r0.get("base_group", r0.get("group", ""))
        c  = base_group_color(bg, base_groups)
        traces.append({
            "type": "scatter", "mode": "lines",
            "x": [c0[0], c1[0]], "y": [c0[1], c1[1]],
            "line": {"color": hex_rgba(c, 0.6), "width": 1.8},
            "showlegend": False, "hoverinfo": "skip",
        })
        traces.append({
            "type": "scatter", "mode": "markers",
            "x": [c0[0]], "y": [c0[1]],
            "marker": {"color": hex_rgba(c, 0.5), "size": 9, "symbol": "circle",
                       "line": {"color": c, "width": 1}},
            "name": f"{p} T0", "showlegend": False,
            "hovertemplate": f"<b>{p} T0</b><extra></extra>",
        })
        traces.append({
            "type": "scatter", "mode": "markers",
            "x": [c1[0]], "y": [c1[1]],
            "marker": {"color": c, "size": 10, "symbol": "circle",
                       "line": {"color": "white", "width": 1}},
            "name": f"{p} T84", "showlegend": False,
            "hovertemplate": f"<b>{p} T84</b><extra></extra>",
        })
    for bg in base_groups:
        c = base_group_color(bg, base_groups)
        traces.append({
            "type": "scatter", "mode": "markers+lines",
            "x": [None], "y": [None],
            "marker": {"color": c, "size": 8},
            "line": {"color": c, "width": 2},
            "name": bg, "showlegend": True,
        })
    return traces, pct1, pct2
