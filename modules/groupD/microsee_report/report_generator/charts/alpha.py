"""charts/alpha.py — alpha-diversity chart builders (strip, box, violin, rarefaction, multi-metric)."""

from __future__ import annotations

from typing import Any

import numpy as np

from .metrics import METRIC_LABELS, metric_value
from .stats_helpers import mannwhitney_p, sig_label, wilcoxon_p
from .utils import group_color, hex_rgba


def build_alpha_strip(
    rows: list[dict[str, Any]],
    groups: list[str],
    metric: str = "shannon",
    taxa: list[str] | None = None,
) -> list[dict[str, Any]]:
    traces: list[dict[str, Any]] = []
    for gi, g in enumerate(groups):
        g_rows = [r for r in rows if r["group"] == g]
        c = group_color(g, groups)
        ys = [metric_value(r, metric, taxa) for r in g_rows]
        xs = [gi + ((i * 1337 + 17) % 100 - 50) / 200 for i in range(len(g_rows))]
        avg = float(np.mean(ys)) if ys else 0.0
        lbl = METRIC_LABELS.get(metric, metric)
        traces.append(
            {
                "type": "scatter",
                "mode": "markers",
                "name": g,
                "x": xs,
                "y": ys,
                "marker": {
                    "color": hex_rgba(c, 0.85),
                    "size": 9,
                    "line": {"width": 1, "color": "white"},
                },
                "text": [r["sample_id"] for r in g_rows],
                "hovertemplate": f"<b>%{{text}}</b><br>{lbl}: %{{y:.3f}}<extra></extra>",
                "showlegend": True,
            }
        )
        traces.append(
            {
                "type": "scatter",
                "mode": "lines",
                "x": [gi - 0.3, gi + 0.3],
                "y": [avg, avg],
                "line": {"color": c, "width": 2.5},
                "showlegend": False,
                "hoverinfo": "skip",
            }
        )
    return traces


def build_alpha_box(
    rows: list[dict[str, Any]],
    groups: list[str],
    metric: str = "shannon",
    taxa: list[str] | None = None,
) -> list[dict[str, Any]]:
    lbl = METRIC_LABELS.get(metric, metric)
    return [
        {
            "type": "box",
            "name": g,
            "y": [metric_value(r, metric, taxa) for r in rows if r["group"] == g],
            "boxpoints": "all",
            "jitter": 0.3,
            "pointpos": -1.8,
            "marker": {"color": hex_rgba(group_color(g, groups), 0.7), "size": 5},
            "line": {"color": group_color(g, groups)},
            "hovertemplate": f"<b>{g}</b><br>{lbl}: %{{y:.3f}}<extra></extra>",
        }
        for g in groups
    ]


def build_alpha_violin(
    rows: list[dict[str, Any]],
    groups: list[str],
    metric: str = "shannon",
    taxa: list[str] | None = None,
) -> list[dict[str, Any]]:
    lbl = METRIC_LABELS.get(metric, metric)
    return [
        {
            "type": "violin",
            "name": g,
            "y": [metric_value(r, metric, taxa) for r in rows if r["group"] == g],
            "box": {"visible": True},
            "meanline": {"visible": True},
            "fillcolor": hex_rgba(group_color(g, groups), 0.4),
            "line": {"color": group_color(g, groups)},
            "hovertemplate": f"<b>{g}</b><br>{lbl}: %{{y:.3f}}<extra></extra>",
        }
        for g in groups
    ]


def _bracket_data(
    rows: list[dict[str, Any]],
    groups: list[str],
    base_groups: list[str],
    metric: str,
    taxa: list[str],
) -> dict[str, Any]:
    """Wilcoxon bracket shapes + annotations for one alpha metric box chart."""
    all_vals = [metric_value(r, metric, taxa) for r in rows]
    if not all_vals:
        return {"shapes": [], "annots": [], "y_max": 1.0}
    y_max = max(all_vals)
    y_span = max(y_max - min(0.0, min(all_vals)), 1e-6)
    gpos = {g: i for i, g in enumerate(groups)}
    shapes: list[dict[str, Any]] = []
    annots: list[dict[str, Any]] = []
    highest = y_max

    def _bracket(x0: int, x1: int, y_b: float, label: str, dashed: bool = False) -> None:
        tick = y_span * 0.035
        col, ls = "#777", ("dot" if dashed else "solid")
        shapes.extend(
            [
                {
                    "type": "line",
                    "x0": x0,
                    "x1": x1,
                    "y0": y_b,
                    "y1": y_b,
                    "xref": "x",
                    "yref": "y",
                    "line": {"color": col, "width": 1, "dash": ls},
                },
                {
                    "type": "line",
                    "x0": x0,
                    "x1": x0,
                    "y0": y_b - tick,
                    "y1": y_b,
                    "xref": "x",
                    "yref": "y",
                    "line": {"color": col, "width": 1},
                },
                {
                    "type": "line",
                    "x0": x1,
                    "x1": x1,
                    "y0": y_b - tick,
                    "y1": y_b,
                    "xref": "x",
                    "yref": "y",
                    "line": {"color": col, "width": 1},
                },
            ]
        )
        annots.append(
            {
                "x": (x0 + x1) / 2,
                "y": y_b + y_span * 0.025,
                "xref": "x",
                "yref": "y",
                "text": label,
                "showarrow": False,
                "font": {"size": 9, "color": "#555"},
            }
        )

    for bg in base_groups:
        t0_g = f"{bg}_T0" if f"{bg}_T0" in gpos else None
        t84_g = f"{bg}_T84" if f"{bg}_T84" in gpos else None
        if not t0_g or not t84_g:
            continue
        patients = sorted(set(r["patient"] for r in rows if r.get("base_group", r["group"]) == bg))
        a: list[float] = []
        b: list[float] = []
        for pat in patients:
            t0v = [
                metric_value(r, metric, taxa)
                for r in rows
                if r["patient"] == pat and (r.get("time") or 0) == 0
            ]
            t84v = [
                metric_value(r, metric, taxa)
                for r in rows
                if r["patient"] == pat and (r.get("time") or 0) > 0
            ]
            if t0v and t84v:
                a.append(t0v[0])
                b.append(t84v[0])
        pval = wilcoxon_p(a, b)
        highest += y_span * 0.14
        _bracket(gpos[t0_g], gpos[t84_g], highest, sig_label(pval))

    if len(base_groups) == 2:
        bg1, bg2 = base_groups
        t84_g1 = f"{bg1}_T84"
        t84_g2 = f"{bg2}_T84"
        if t84_g1 in gpos and t84_g2 in gpos:
            v1 = [
                metric_value(r, metric, taxa)
                for r in rows
                if r.get("base_group", r["group"]) == bg1 and r.get("timepoint") == "T84"
            ]
            v2 = [
                metric_value(r, metric, taxa)
                for r in rows
                if r.get("base_group", r["group"]) == bg2 and r.get("timepoint") == "T84"
            ]
            p = mannwhitney_p(v1, v2)
            highest += y_span * 0.16
            _bracket(gpos[t84_g1], gpos[t84_g2], highest, sig_label(p) + " (T84)", dashed=True)

    return {"shapes": shapes, "annots": annots, "y_max": round(highest + y_span * 0.10, 4)}


def build_all_alpha_metrics(
    rows: list[dict[str, Any]],
    groups: list[str],
    taxa: list[str],
    base_groups: list[str] | None = None,
) -> dict[str, Any]:
    """Pre-compute strip/box/violin + significance brackets for every metric."""
    bgs = base_groups or sorted(set(r.get("base_group", r["group"]) for r in rows))
    return {
        metric: {
            "strip": build_alpha_strip(rows, groups, metric, taxa),
            "box": build_alpha_box(rows, groups, metric, taxa),
            "violin": build_alpha_violin(rows, groups, metric, taxa),
            "brackets": _bracket_data(rows, groups, bgs, metric, taxa),
        }
        for metric in METRIC_LABELS
    }


def build_rarefaction(
    rows: list[dict[str, Any]], taxa: list[str], groups: list[str]
) -> list[dict[str, Any]]:
    """Expected species accumulation ± 1 SD shaded band per group."""
    depths = [10, 25, 50, 100, 200, 500, 1_000, 2_000, 5_000, 10_000]
    traces: list[dict[str, Any]] = []
    for g in groups:
        g_rows = [r for r in rows if r["group"] == g]
        c = group_color(g, groups)
        probs_list = [np.array([float(r.get(t) or 0) / 100.0 for t in taxa]) for r in g_rows]
        probs_list = [p / p.sum() if p.sum() > 0 else p for p in probs_list]
        curves = [[float(np.sum(1 - (1 - p) ** d)) for d in depths] for p in probs_list]
        if not curves:
            continue
        mean_c = np.mean(curves, axis=0)
        std_c = np.std(curves, axis=0) if len(curves) > 1 else np.zeros(len(depths))
        upper = [round(float(m + s), 2) for m, s in zip(mean_c, std_c, strict=False)]
        lower = [round(max(0.0, float(m - s)), 2) for m, s in zip(mean_c, std_c, strict=False)]
        mean_v = [round(float(v), 2) for v in mean_c]
        sd_avg = round(float(np.mean(std_c)), 1)
        # Upper boundary (invisible — tonexty fills from below)
        traces.append(
            {
                "type": "scatter",
                "mode": "lines",
                "x": depths,
                "y": upper,
                "line": {"width": 0},
                "showlegend": False,
                "hoverinfo": "skip",
                "name": g,
            }
        )
        # Lower boundary fills up to upper
        traces.append(
            {
                "type": "scatter",
                "mode": "lines",
                "x": depths,
                "y": lower,
                "fill": "tonexty",
                "fillcolor": hex_rgba(c, 0.15),
                "line": {"width": 0},
                "showlegend": False,
                "hoverinfo": "skip",
                "name": g,
            }
        )
        # Mean line
        traces.append(
            {
                "type": "scatter",
                "mode": "lines",
                "name": g,
                "x": depths,
                "y": mean_v,
                "line": {"color": c, "width": 2.5},
                "hovertemplate": (
                    f"<b>{g}</b><br>Depth: %{{x}}<br>"
                    f"Expected taxa: %{{y:.1f}} ± {sd_avg}<extra></extra>"
                ),
            }
        )
    return traces


def build_multimet_alpha(
    rows: list[dict[str, Any]], taxa: list[str], groups: list[str]
) -> list[dict[str, Any]]:
    traces: list[dict[str, Any]] = []
    for g in groups:
        g_rows = [r for r in rows if r["group"] == g]
        c = group_color(g, groups)
        traces.append(
            {
                "type": "bar",
                "name": f"{g} — Observed",
                "x": [r["sample_id"] for r in g_rows],
                "y": [metric_value(r, "observed", taxa) for r in g_rows],
                "marker": {"color": hex_rgba(c, 0.8)},
                "hovertemplate": f"<b>{g}</b><br>Observed taxa: %{{y}}<extra></extra>",
                "yaxis": "y",
            }
        )
        traces.append(
            {
                "type": "scatter",
                "mode": "markers",
                "name": f"{g} — Pielou J′",
                "x": [r["sample_id"] for r in g_rows],
                "y": [metric_value(r, "pielou", taxa) for r in g_rows],
                "marker": {"color": c, "size": 8, "symbol": "diamond"},
                "hovertemplate": f"<b>{g}</b><br>Pielou J′: %{{y:.3f}}<extra></extra>",
                "yaxis": "y2",
            }
        )
    return traces
