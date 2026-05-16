"""charts/clinical.py — clinical outcome chart builders (slopegraph, Shannon correlation)."""

from __future__ import annotations
import math
import numpy as np
from .utils import hex_rgba, _base_group_color, _sorted_timepoints


def _pearson_p(r_val: float, n: int) -> float:
    """Two-tailed p-value for Pearson r using Fisher z-transformation."""
    if n < 4:
        return 1.0
    z   = math.atanh(min(max(r_val, -0.9999), 0.9999)) * math.sqrt(n - 3)
    az  = abs(z)
    t   = 1.0 / (1.0 + 0.2316419 * az)
    pd  = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))))
    p1  = pd * math.exp(-0.5 * az ** 2) / math.sqrt(2 * math.pi)
    return min(1.0, 2 * p1)


def _spearman_r(a: list[float], b: list[float]) -> tuple[float, float]:
    """Spearman rho and two-tailed p-value."""
    n = len(a)
    if n < 3:
        return 0.0, 1.0
    def _ranks(v: list[float]) -> list[float]:
        sv = sorted(enumerate(v), key=lambda x: x[1])
        r  = [0.0] * n
        i  = 0
        while i < n:
            j = i
            while j < n - 1 and sv[j][1] == sv[j + 1][1]:
                j += 1
            avg = (i + j + 2) / 2.0
            for k in range(i, j + 1):
                r[sv[k][0]] = avg
            i = j + 1
        return r
    ra, rb = _ranks(a), _ranks(b)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    da  = math.sqrt(sum((x - ma) ** 2 for x in ra))
    db  = math.sqrt(sum((y - mb) ** 2 for y in rb))
    rho = round(num / (da * db), 3) if da and db else 0.0
    p   = _pearson_p(rho, n)
    return rho, p


def build_taxa_clinical_heatmap(rows: list[dict], taxa: list[str]) -> dict:
    """Spearman rho heatmap: all taxa × all clinical variables.

    Uses the change (Δ = T84 − T0) for both taxon abundances and clinical values,
    matching the paper's approach of asking 'do microbiome shifts correlate with
    clinical improvements?'
    """
    clinical_fields = [f for f in ("sixmwt", "il18") if any(float(r.get(f) or 0) != 0 for r in rows)]
    if not clinical_fields:
        return {"z": [], "x": [], "y": [], "text": [], "p": []}

    patients = sorted(set(r["patient"] for r in rows))

    # Compute Δ per patient for each taxon and clinical field
    delta_tax: dict[str, list[float]] = {t: [] for t in taxa}
    delta_clin: dict[str, list[float]] = {f: [] for f in clinical_fields}
    valid_patients = []

    for p in patients:
        t0r  = [r for r in rows if r["patient"] == p and (r.get("time") or 0) == 0]
        t84r = sorted([r for r in rows if r["patient"] == p and (r.get("time") or 0) > 0],
                      key=lambda r: r.get("time") or 0, reverse=True)
        if not t0r or not t84r:
            continue
        r0, r84 = t0r[0], t84r[0]
        valid_patients.append(p)
        for t in taxa:
            delta_tax[t].append(float(r84.get(t) or 0) - float(r0.get(t) or 0))
        for f in clinical_fields:
            delta_clin[f].append(float(r84.get(f) or 0) - float(r0.get(f) or 0))

    if len(valid_patients) < 3:
        return {"z": [], "x": [], "y": [], "text": [], "p": []}

    field_labels = {"sixmwt": "Δ 6MWT (m)", "il18": "Δ IL-18 (pg/mL)"}
    x_labels = [field_labels.get(f, f) for f in clinical_fields]
    z, text_mat, p_mat = [], [], []

    for t in taxa:
        row_z, row_t, row_p = [], [], []
        for f in clinical_fields:
            rho, p = _spearman_r(delta_tax[t], delta_clin[f])
            row_z.append(rho)
            row_p.append(round(p, 3))
            stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            row_t.append(f"ρ={rho:+.2f}{stars}")
        z.append(row_z); text_mat.append(row_t); p_mat.append(row_p)

    return {
        "z": z, "x": x_labels, "y": list(taxa),
        "text": text_mat, "p": p_mat,
    }


def build_clinical_slope(rows: list[dict], field: str) -> list[dict]:
    patients    = sorted(set(r["patient"] for r in rows))
    timepoints  = _sorted_timepoints(rows)
    base_groups = sorted(set(r.get("base_group", r["group"]) for r in rows))
    traces: list[dict] = []

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
            "y": [float(r.get(field) or 0) for r in p_rows],
            "line": {"color": hex_rgba(c, 0.55), "width": 1.5},
            "marker": {"size": 6, "color": c},
            "showlegend": False,
            "hovertemplate": f"<b>{p}</b><br>%{{x}}<br>%{{y:.1f}}<extra></extra>",
        })

    for bg in base_groups:
        c = _base_group_color(bg, base_groups)
        means, xs_used = [], []
        for tp in timepoints:
            vals = [float(r.get(field) or 0)
                    for r in rows
                    if r.get("base_group", r["group"]) == bg and r["timepoint"] == tp]
            if vals:
                means.append(round(float(np.mean(vals)), 2))
                xs_used.append(tp)
        if means:
            traces.append({
                "type": "scatter", "mode": "lines+markers",
                "name": f"{bg} mean", "x": xs_used, "y": means,
                "line": {"color": c, "width": 3, "dash": "dash"},
                "marker": {"size": 10, "color": c, "symbol": "diamond"},
                "hovertemplate": f"<b>{bg} mean</b><br>%{{x}}<br>%{{y:.1f}}<extra></extra>",
            })
    return traces


def build_clinical_correlation(rows: list[dict], field: str) -> tuple[list[dict], float, float]:
    """Shannon vs clinical field scatter + regression line. Returns (traces, r, p)."""
    base_groups = sorted(set(r.get("base_group", r["group"]) for r in rows))
    valid = [r for r in rows if float(r.get(field) or 0) > 0 and float(r.get("shannon") or 0) > 0]
    if len(valid) < 3:
        return [], 0.0, 1.0

    xs = np.array([float(r.get("shannon") or 0) for r in valid])
    ys = np.array([float(r.get(field) or 0) for r in valid])
    mx, my = float(np.mean(xs)), float(np.mean(ys))
    num   = float(np.sum((xs - mx) * (ys - my)))
    dx    = float(np.sqrt(np.sum((xs - mx) ** 2)))
    dy    = float(np.sqrt(np.sum((ys - my) ** 2)))
    r_val = round(num / (dx * dy), 3) if dx and dy else 0.0
    p_val = round(_pearson_p(r_val, len(valid)), 3)

    slope     = num / (dx ** 2 + 1e-10) if dx else 0.0
    intercept = my - slope * mx
    x_range   = [float(np.min(xs)), float(np.max(xs))]

    traces: list[dict] = []
    for bg in base_groups:
        c   = _base_group_color(bg, base_groups)
        pts = [r for r in valid if r.get("base_group", r.get("group")) == bg]
        if not pts:
            continue
        traces.append({
            "type": "scatter", "mode": "markers", "name": bg,
            "x": [float(r.get("shannon") or 0) for r in pts],
            "y": [float(r.get(field) or 0) for r in pts],
            "text": [r["sample_id"] for r in pts],
            "marker": {"color": hex_rgba(c, 0.85), "size": 9,
                       "line": {"width": 1, "color": "white"}},
            "hovertemplate": "<b>%{text}</b><br>%{x:.3f} / %{y:.1f}<extra></extra>",
        })
    traces.append({
        "type": "scatter", "mode": "lines", "name": "Regression",
        "x": x_range,
        "y": [slope * x + intercept for x in x_range],
        "line": {"color": hex_rgba("#C4A0A8", 0.7), "width": 1.5, "dash": "dash"},
        "showlegend": False, "hoverinfo": "skip",
    })
    return traces, r_val, p_val
