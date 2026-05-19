"""charts/taxonomy.py — stacked bar, top-N, donut, and sunburst chart builders."""

from __future__ import annotations

from typing import Any

import numpy as np

from .utils import group_color, hex_rgba, taxon_color


def build_taxonomy_stacked(
    rows: list[dict[str, Any]],
    taxa: list[str],
) -> list[dict[str, Any]]:
    """Stacked bar chart: one trace per taxon, x = sample IDs."""
    sample_ids = [r["sample_id"] for r in rows]
    return [
        {
            "type": "bar",
            "name": taxon,
            "x": sample_ids,
            "y": [round(float(r.get(taxon) or 0), 2) for r in rows],
            "marker": {"color": hex_rgba(taxon_color(taxon), 0.9)},
            "hovertemplate": f"<b>{taxon}</b><br>%{{x}}<br>%{{y:.1f}}%<extra></extra>",
        }
        for taxon in taxa
    ]


def build_top_taxa(rows: list[dict[str, Any]], taxa: list[str]) -> list[dict[str, Any]]:
    """Horizontal bar chart of taxa sorted by cohort-mean relative abundance."""
    means = {t: float(np.mean([float(r.get(t) or 0) for r in rows])) for t in taxa}
    ordered = sorted(means, key=lambda t: means[t], reverse=True)
    return [
        {
            "type": "bar",
            "orientation": "h",
            "x": [round(means[t], 2) for t in ordered],
            "y": ordered,
            "marker": {"color": [taxon_color(t) for t in ordered]},
            "hovertemplate": "<b>%{y}</b><br>%{x:.2f}%<extra></extra>",
        }
    ]


def build_donut(
    rows: list[dict[str, Any]],
    taxa: list[str],
    groups: list[str],
) -> list[dict[str, Any]]:
    """One donut per group showing mean relative taxon abundance."""
    col_w = 1.0 / len(groups)
    traces: list[dict[str, Any]] = []
    for gi, g in enumerate(groups):
        g_rows = [r for r in rows if r["group"] == g]
        if not g_rows:
            continue
        traces.append(
            {
                "type": "pie",
                "hole": 0.45,
                "name": g,
                "labels": list(taxa),
                "values": [float(np.mean([float(r.get(t) or 0) for r in g_rows])) for t in taxa],
                "marker": {"colors": [taxon_color(t) for t in taxa]},
                "domain": {"x": [gi * col_w, (gi + 1) * col_w - 0.02], "y": [0.1, 1.0]},
                "title": {"text": g, "font": {"size": 11}},
                "hovertemplate": "<b>%{label}</b><br>%{percent:.1%}<extra></extra>",
            }
        )
    return traces


def build_taxonomy_views(
    rows: list[dict[str, Any]],
    taxa: list[str],
    base_groups: list[str],
) -> dict[str, Any]:
    """Pre-compute stacked-bar traces for every tp x group x topn combination."""
    views: dict[str, list[dict[str, Any]]] = {}
    tp_map: dict[str, str | None] = {"both": None, "T0": "T0", "T84": "T84"}
    grp_map: dict[str, str | None] = {"all": None, **{bg: bg for bg in base_groups}}
    topns: list[int | None] = [5, 10, None]

    for tp_key, tp_val in tp_map.items():
        for grp_key, grp_val in grp_map.items():
            filtered = rows
            if tp_val:
                filtered = [r for r in filtered if r.get("timepoint") == tp_val]
            if grp_val:
                filtered = [r for r in filtered if r.get("base_group", r.get("group")) == grp_val]
            if not filtered:
                continue
            for topn in topns:
                topn_key = "all" if topn is None else str(topn)
                top_taxa = taxa[:topn] if topn else taxa
                views[f"{tp_key}:{grp_key}:{topn_key}"] = build_taxonomy_stacked(
                    filtered,
                    top_taxa,
                )

    return views


def build_sunburst(
    rows: list[dict[str, Any]],
    taxa: list[str],
    groups: list[str],
) -> list[dict[str, Any]]:
    """Sunburst chart: group as root, taxa as leaves weighted by mean abundance."""
    labels: list[str] = []
    parents: list[str] = []
    values: list[float] = []
    colors: list[str] = []

    for g in groups:
        g_rows = [r for r in rows if r["group"] == g]
        child_vals = [
            round(float(np.mean([float(r.get(t) or 0) for r in g_rows])), 1) for t in taxa
        ]
        g_total = round(sum(child_vals), 1)
        labels.append(g)
        parents.append("")
        values.append(g_total)
        colors.append(group_color(g, groups))
        for t, v in zip(taxa, child_vals, strict=False):
            if v > 0:
                labels.append(f"{g}:{t}")
                parents.append(g)
                values.append(v)
                colors.append(taxon_color(t))

    return [
        {
            "type": "sunburst",
            "labels": labels,
            "parents": parents,
            "values": values,
            "marker": {"colors": colors},
            "branchvalues": "total",
            "hovertemplate": "<b>%{label}</b><br>%{value:.1f}%<extra></extra>",
            "textinfo": "label+percent entry",
        }
    ]
