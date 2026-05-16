"""charts/individual.py — per-patient chart builders (slopegraph, stability, rank, radar, faceted, NMDS trajectories)."""

from __future__ import annotations
import numpy as np
from .utils import hex_rgba, taxon_color, _base_group_color, _sorted_timepoints
from .distances import rows_to_ab, bray_curtis_matrix, pcoa


def build_paired_slope(rows: list[dict], groups: list[str]) -> list[dict]:
    """Per-patient Shannon T0→T84 slopegraph."""
    timepoints  = _sorted_timepoints(rows)
    base_groups = sorted(set(r.get("base_group", r["group"]) for r in rows))
    traces: list[dict] = []
    patients = sorted(set(r["patient"] for r in rows))
    for p in patients:
        p_rows = sorted([r for r in rows if r["patient"] == p],
                        key=lambda r: r.get("time") or 0)
        if len(p_rows) < 2:
            continue
        bg = p_rows[0].get("base_group", p_rows[0]["group"])
        c  = _base_group_color(bg, base_groups)
        traces.append({
            "type": "scatter", "mode": "lines+markers", "name": p,
            "x": [r["timepoint"] for r in p_rows],
            "y": [float(r["shannon"]) for r in p_rows],
            "line": {"color": hex_rgba(c, 0.55), "width": 1.5},
            "marker": {"size": 7, "color": c},
            "showlegend": False,
            "hovertemplate": f"<b>{p}</b><br>%{{x}}<br>Shannon: %{{y:.3f}}<extra></extra>",
        })
    for bg in base_groups:
        c = _base_group_color(bg, base_groups)
        means, xs_used = [], []
        for tp in timepoints:
            vals = [float(r["shannon"]) for r in rows
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
                "hovertemplate": f"<b>{bg} mean</b><br>%{{x}}<br>Shannon: %{{y:.3f}}<extra></extra>",
            })
    return traces


def build_stability_bar(rows: list[dict], taxa: list[str]) -> list[dict]:
    """Per-patient Bray-Curtis dissimilarity T0 vs T84, sorted horizontal bar."""
    base_groups = sorted(set(r.get("base_group", r["group"]) for r in rows))
    patients    = sorted(set(r["patient"] for r in rows))
    scores = []
    for p in patients:
        r0_list  = [r for r in rows if r["patient"] == p and (r.get("time") or 0) == 0]
        r84_list = sorted([r for r in rows if r["patient"] == p and (r.get("time") or 0) > 0],
                          key=lambda r: r.get("time") or 0, reverse=True)
        if not r0_list or not r84_list:
            continue
        r0, r84 = r0_list[0], r84_list[0]
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
        "marker": {"color": [hex_rgba(_base_group_color(s["base_group"], base_groups), 0.8) for s in scores]},
        "text": [str(s["bc"]) for s in scores],
        "textposition": "outside",
        "hovertemplate": "<b>%{y}</b><br>BC dissimilarity: %{x:.3f}<extra></extra>",
    }]


def build_diversity_rank(rows: list[dict]) -> list[dict]:
    """Samples ranked low→high by Shannon; circle = T0, diamond = T84."""
    base_groups = sorted(set(r.get("base_group", r["group"]) for r in rows))
    sorted_rows = sorted(rows, key=lambda r: float(r.get("shannon") or 0))
    traces = []
    for bg in base_groups:
        c    = _base_group_color(bg, base_groups)
        pts  = [(i, r) for i, r in enumerate(sorted_rows)
                if r.get("base_group", r.get("group")) == bg]
        t0   = [(i, r) for i, r in pts if (r.get("time") or 0) == 0]
        t84  = [(i, r) for i, r in pts if (r.get("time") or 0) > 0]
        if t0:
            traces.append({
                "type": "scatter", "mode": "markers", "name": f"{bg} T0",
                "x": [i for i, _ in t0],
                "y": [float(r.get("shannon") or 0) for _, r in t0],
                "text": [r["sample_id"] for _, r in t0],
                "marker": {"color": hex_rgba(c, 0.6), "size": 8, "symbol": "circle",
                           "line": {"width": 1, "color": "white"}},
                "hovertemplate": "<b>%{text}</b>: %{y:.3f}<extra></extra>",
            })
        if t84:
            traces.append({
                "type": "scatter", "mode": "markers", "name": f"{bg} T84",
                "x": [i for i, _ in t84],
                "y": [float(r.get("shannon") or 0) for _, r in t84],
                "text": [r["sample_id"] for _, r in t84],
                "marker": {"color": c, "size": 9, "symbol": "diamond",
                           "line": {"width": 1, "color": "white"}},
                "hovertemplate": "<b>%{text}</b>: %{y:.3f}<extra></extra>",
            })
    return traces


def build_patient_radar(rows: list[dict], taxa: list[str]) -> list[dict]:
    """Radar: group mean T0 (filled) vs group mean T84 (dashed) per base_group."""
    base_groups = sorted(set(r.get("base_group", r["group"]) for r in rows))
    short_taxa  = [t.replace("aceae", "")[:10] for t in taxa]
    traces = []
    for bg in base_groups:
        c        = _base_group_color(bg, base_groups)
        t0_rows  = [r for r in rows if r.get("base_group", r.get("group")) == bg and (r.get("time") or 0) == 0]
        t84_rows = [r for r in rows if r.get("base_group", r.get("group")) == bg and (r.get("time") or 0) > 0]
        for label, g_rows, filled in [("T0", t0_rows, True), ("T84", t84_rows, False)]:
            if not g_rows:
                continue
            vals = [round(float(np.mean([float(r.get(t) or 0) for r in g_rows])), 1) for t in taxa]
            traces.append({
                "type": "scatterpolar",
                "fill": "toself" if filled else "none",
                "r": vals + [vals[0]],
                "theta": short_taxa + [short_taxa[0]],
                "fillcolor": hex_rgba(c, 0.2) if filled else "rgba(0,0,0,0)",
                "line": {"color": hex_rgba(c, 0.8) if filled else c,
                         "width": 2 if filled else 2.5,
                         "dash": "solid" if filled else "dash"},
                "name": f"{bg} {label}",
                "hovertemplate": "<b>%{theta}</b>: %{r:.1f}%<extra></extra>",
            })
    return traces


def build_faceted_composition(rows: list[dict], taxa: list[str]) -> dict:
    """True small-multiples: one subplot per patient, T0 vs T84 stacked bars."""
    import math
    patients = sorted(set(r["patient"] for r in rows))
    n        = len(patients)
    n_cols   = min(4, n)
    n_rows   = math.ceil(n / n_cols)
    col_gap, row_gap = 0.04, 0.10
    col_w = (1.0 - col_gap * (n_cols - 1)) / n_cols
    row_h = (1.0 - row_gap * (n_rows - 1)) / n_rows

    data: list[dict] = []
    annotations: list[dict] = []
    layout: dict = {"barmode": "stack", "showlegend": True,
                    "height": max(360, n_rows * 240),
                    "margin": {"l": 40, "r": 10, "t": 20, "b": 40}}
    shown: set = set()

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

        r0_list  = [r for r in rows if r["patient"] == p and (r.get("time") or 0) == 0]
        r84_list = sorted([r for r in rows if r["patient"] == p and (r.get("time") or 0) > 0],
                          key=lambda r: r.get("time") or 0, reverse=True)

        for r, tp in [(r0_list[0] if r0_list else None, "T0"),
                      (r84_list[0] if r84_list else None, "T84")]:
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


def build_nmds_trajectories(rows: list[dict], taxa: list[str]) -> tuple[list[dict], float, float]:
    """PCoA (Bray-Curtis) with T0→T84 arrows per patient."""
    base_groups = sorted(set(r.get("base_group", r["group"]) for r in rows))
    ab  = rows_to_ab(rows, taxa)
    mat = bray_curtis_matrix(ab)
    xs, ys, pct1, pct2 = pcoa(mat)
    coord_map = {r["sample_id"]: (xs[i], ys[i]) for i, r in enumerate(rows)}
    patients  = sorted(set(r["patient"] for r in rows))
    traces: list[dict] = []
    for p in patients:
        r0_list  = [r for r in rows if r["patient"] == p and (r.get("time") or 0) == 0]
        r84_list = sorted([r for r in rows if r["patient"] == p and (r.get("time") or 0) > 0],
                          key=lambda r: r.get("time") or 0, reverse=True)
        if not r0_list or not r84_list:
            continue
        r0, r84 = r0_list[0], r84_list[0]
        c0 = coord_map.get(r0["sample_id"])
        c1 = coord_map.get(r84["sample_id"])
        if not c0 or not c1:
            continue
        bg = r0.get("base_group", r0.get("group", ""))
        c  = _base_group_color(bg, base_groups)
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
        c = _base_group_color(bg, base_groups)
        traces.append({
            "type": "scatter", "mode": "markers+lines",
            "x": [None], "y": [None],
            "marker": {"color": c, "size": 8},
            "line": {"color": c, "width": 2},
            "name": bg, "showlegend": True,
        })
    return traces, pct1, pct2
