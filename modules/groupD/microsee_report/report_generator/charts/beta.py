"""charts/beta.py — beta-diversity chart builders.

PCoA, NMDS, dendrogram, and delta abundance heatmap.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .distances import (
    average_linkage_dendrogram,
    bray_curtis_matrix,
    jaccard_matrix,
    pcoa,
    rows_to_ab,
)
from .preprocessing import get_patient_timepoints, get_unique_patients
from .utils import group_color, hex_rgba

_MIN_DENDROGRAM_SAMPLES = 3


def build_pcoa_chart(
    rows: list[dict[str, Any]],
    taxa: list[str],
    groups: list[str],
    dist_type: str,
    bc_matrix: np.ndarray | None = None,
) -> tuple[list[dict[str, Any]], float, float]:
    """PCoA scatter coloured by group (Bray-Curtis or Jaccard distance)."""
    if dist_type == "bray" and bc_matrix is not None:
        mat = bc_matrix
    else:
        ab = rows_to_ab(rows, taxa)
        mat = bray_curtis_matrix(ab) if dist_type == "bray" else jaccard_matrix(ab)
    xs, ys, pct1, pct2 = pcoa(mat)
    traces: list[dict[str, Any]] = []
    for g in groups:
        idxs = [i for i, r in enumerate(rows) if r["group"] == g]
        c = group_color(g, groups)
        traces.append(
            {
                "type": "scatter",
                "mode": "markers",
                "name": g,
                "x": [xs[i] for i in idxs],
                "y": [ys[i] for i in idxs],
                "text": [rows[i]["sample_id"] for i in idxs],
                "marker": {
                    "color": hex_rgba(c, 0.85),
                    "size": 10,
                    "line": {"width": 1, "color": "white"},
                },
                "hovertemplate": "<b>%{text}</b><extra></extra>",
            }
        )
    return traces, pct1, pct2


def build_nmds_plot(
    rows: list[dict[str, Any]],
    taxa: list[str],
    groups: list[str],
    bc_matrix: np.ndarray | None = None,
) -> tuple[list[dict[str, Any]], float, float]:
    """NMDS-style plot (PCoA Bray-Curtis, diamond markers — matches NMDSPlot.tsx)."""
    mat = bc_matrix if bc_matrix is not None else bray_curtis_matrix(rows_to_ab(rows, taxa))
    xs, ys, pct1, pct2 = pcoa(mat)
    traces: list[dict[str, Any]] = []
    for g in groups:
        idxs = [i for i, r in enumerate(rows) if r["group"] == g]
        c = group_color(g, groups)
        traces.append(
            {
                "type": "scatter",
                "mode": "markers",
                "name": g,
                "x": [xs[i] for i in idxs],
                "y": [ys[i] for i in idxs],
                "text": [rows[i]["sample_id"] for i in idxs],
                "marker": {
                    "color": hex_rgba(c, 0.85),
                    "size": 10,
                    "symbol": "diamond",
                    "line": {"width": 1, "color": "white"},
                },
                "hovertemplate": "<b>%{text}</b><extra></extra>",
            }
        )
    return traces, pct1, pct2


def build_dendrogram(
    rows: list[dict[str, Any]],
    taxa: list[str],
    groups: list[str],
    bc_matrix: np.ndarray | None = None,
) -> list[dict[str, Any]]:
    """Hierarchical clustering dendrogram (average-linkage, Bray-Curtis)."""
    n = len(rows)
    if n < _MIN_DENDROGRAM_SAMPLES:
        return []
    mat = bc_matrix if bc_matrix is not None else bray_curtis_matrix(rows_to_ab(rows, taxa))
    labels = [r["sample_id"] for r in rows]
    sid_grp = {r["sample_id"]: r.get("group", "") for r in rows}
    segments = average_linkage_dendrogram(mat)

    traces: list[dict[str, Any]] = [
        {
            "type": "scatter",
            "mode": "lines",
            "x": xs,
            "y": ys,
            "line": {"color": "#C4A0A8", "width": 1.5},
            "showlegend": False,
            "hoverinfo": "skip",
        }
        for xs, ys in segments
    ]

    traces.append(
        {
            "type": "scatter",
            "mode": "markers",
            "x": [0.0] * n,
            "y": list(range(n)),
            "marker": {
                "color": [group_color(sid_grp.get(sid, ""), groups) for sid in labels],
                "size": 8,
            },
            "text": labels,
            "showlegend": False,
            "hovertemplate": "<b>%{text}</b><extra></extra>",
        }
    )
    traces += [
        {
            "type": "scatter",
            "mode": "markers",
            "x": [None],
            "y": [None],
            "marker": {"color": group_color(g, groups), "size": 8},
            "name": g,
            "showlegend": True,
        }
        for g in groups
    ]
    return traces


def build_delta_heatmap(rows: list[dict[str, Any]], taxa: list[str]) -> list[dict[str, Any]]:
    """Delta abundance heatmap: (T84 - T0) normalised to % per patient per taxon."""
    patients = get_unique_patients(rows)
    patient_labels: list[str] = []
    delta_rows: list[list[float]] = []

    for p in patients:
        r0, r84 = get_patient_timepoints(rows, p)
        if r0 is None or r84 is None:
            continue
        tot0 = sum(float(r0.get(t) or 0) for t in taxa) or 1
        tot84 = sum(float(r84.get(t) or 0) for t in taxa) or 1
        delta_rows.append(
            [
                round(
                    float(r84.get(t) or 0) / tot84 * 100 - float(r0.get(t) or 0) / tot0 * 100,
                    1,
                )
                for t in taxa
            ]
        )
        patient_labels.append(p)

    if not delta_rows:
        return []

    arr = np.array(delta_rows)
    max_abs = np.max(np.abs(arr), axis=0)
    order: list[int] = [int(i) for i in np.argsort(max_abs)[::-1]]
    sorted_taxa = [taxa[i] for i in order]
    z = [[delta_rows[pi][i] for pi in range(len(patient_labels))] for i in order]

    return [
        {
            "type": "heatmap",
            "z": z,
            "x": patient_labels,
            "y": sorted_taxa,
            "colorscale": [
                [0.0, "#2050A0"],
                [0.35, "#A0C0F8"],
                [0.5, "#FFF8F4"],
                [0.65, "#F8C0A0"],
                [1.0, "#C03030"],
            ],
            "zmid": 0,
            "zmin": -12,
            "zmax": 12,
            "colorbar": {"title": {"text": "Δ%"}, "len": 0.8, "thickness": 14},
            "xgap": 1.5,
            "ygap": 1.5,
            "hovertemplate": "<b>%{x}</b> · <b>%{y}</b><br>Δ = %{z:+.1f}%<extra></extra>",
        }
    ]
