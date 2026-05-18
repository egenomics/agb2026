"""charts/stats.py — statistical summary builders (longitudinal, diversity table, PERMANOVA,
                       Wilcoxon / Mann-Whitney tests, LME-style trajectory)."""
from __future__ import annotations
from typing import Any
import hashlib
import math
import numpy as np
from .utils import base_group_color, hex_rgba
from .preprocessing import sorted_timepoints, get_base_groups
from .metrics import METRIC_LABELS, metric_value
from .stats_helpers import wilcoxon_p, sig_label
from .distances import rows_to_ab, bray_curtis_matrix


def build_longitudinal(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Group mean Shannon per timepoint, connected across time."""
    timepoints  = sorted_timepoints(rows)
    base_groups = get_base_groups(rows)
    traces = []
    for bg in base_groups:
        c = base_group_color(bg, base_groups)
        xs_used, ys = [], []
        for tp in timepoints:
            vals = [float(r["shannon"])
                    for r in rows
                    if r.get("base_group", r["group"]) == bg and r["timepoint"] == tp]
            if vals:
                ys.append(round(float(np.mean(vals)), 4))
                xs_used.append(tp)
        if ys:
            traces.append({
                "type": "scatter", "mode": "lines+markers", "name": bg,
                "x": xs_used, "y": ys,
                "line": {"color": c, "width": 2},
                "marker": {"size": 8, "color": c},
                "hovertemplate": f"<b>{bg}</b><br>%{{x}}<br>Shannon: %{{y:.3f}}<extra></extra>",
            })
    return traces


def build_stats_table(rows: list[dict[str, Any]], groups: list[str], taxa: list[str] | None = None) -> dict[str, Any]:
    metrics = list(METRIC_LABELS.keys())
    header = ["Group", "N"] + [f"{METRIC_LABELS[m]} mean ± SD" for m in metrics]
    table_rows = []
    for g in groups:
        g_rows = [r for r in rows if r["group"] == g]
        cells: list = [g, len(g_rows)]
        for m in metrics:
            vals = [metric_value(r, m, taxa) for r in g_rows]
            non_zero = [v for v in vals if v > 0]
            if non_zero:
                cells.append(f"{np.mean(non_zero):.3f} ± {np.std(non_zero):.3f}")
            else:
                cells.append("—")
        table_rows.append(cells)
    return {"header": header, "rows": table_rows}


def build_permanova_table(
    rows: list[dict[str, Any]],
    taxa: list[str],
    bc_matrix: "np.ndarray | None" = None,
) -> dict[str, Any]:
    """PERMANOVA: pseudo-F, R², p-value (99 permutations, Bray-Curtis)."""
    n = len(rows)
    if n < 4 or not taxa:
        return {"header": [], "rows": [], "top_name": "", "top_R2": 0.0, "top_is_individual": False}

    if bc_matrix is None:
        bc_matrix = bray_curtis_matrix(rows_to_ab(rows, taxa))
    mat = bc_matrix

    # Pre-compute upper-triangle squared distances — invariant under label permutations
    ui, uj = np.triu_indices(n, k=1)
    d2_all = mat[ui, uj] ** 2
    ss_tot = float(d2_all.sum()) / n

    def pseudo_f(labels_arr: np.ndarray) -> tuple[float, float]:
        gs   = np.unique(labels_arr)
        ss_w = 0.0
        for g in gs:
            idx = np.where(labels_arr == g)[0]
            ng  = len(idx)
            if ng < 2:
                continue
            sub        = mat[np.ix_(idx, idx)]
            si, sj     = np.triu_indices(ng, k=1)
            ss_w      += float((sub[si, sj] ** 2).sum()) / ng
        ss_b = ss_tot - ss_w
        k    = len(gs)
        R2   = ss_b / (ss_tot or 1.0)
        F    = (ss_b / max(k - 1, 1)) / max(ss_w / max(n - k, 1), 1e-10)
        return F, R2

    _key = ",".join(sorted(r.get("sample_id", "") for r in rows)).encode()
    seed = int(hashlib.md5(_key).hexdigest()[:8], 16) or 42
    rng  = np.random.default_rng(seed)

    def permute(labels: list[str], n_perm: int = 99) -> tuple[float, float, float]:
        labels_arr = np.array(labels)
        obs_F, obs_R2 = pseudo_f(labels_arr)
        cnt = 1
        for _ in range(n_perm):
            perm = rng.permutation(n)
            if pseudo_f(labels_arr[perm])[0] >= obs_F:
                cnt += 1
        return obs_F, obs_R2, cnt / (n_perm + 1)

    def _bg(r: dict) -> str:
        return str(r.get("base_group") or r.get("group") or "").replace("_T0", "").replace("_T84", "")
    def _pid(r: dict) -> str:
        return str(r.get("patient") or r.get("sample_id") or "")

    tests = [
        ("Supplementation Group", [_bg(r) for r in rows]),
        ("Timepoint",             [str(r.get("time") or r.get("timepoint") or "?") for r in rows]),
        ("Full Subgroup",         [str(r.get("group") or "") for r in rows]),
        ("Individual",            [_pid(r) for r in rows]),
    ]

    results = sorted([
        {"name": name, "R2": round(R2, 3), "F": round(F, 2), "p": round(p, 3)}
        for name, labels in tests
        for F, R2, p in [permute(labels)]
    ], key=lambda x: x["R2"], reverse=True)
    top = results[0]

    return {
        "header": ["Variable", "R²", "Pseudo-F", "p-value", "Result"],
        "rows": [
            [r["name"], r["R2"], r["F"], r["p"], "✓ Significant" if r["p"] < 0.05 else "Not significant"]
            for r in results
        ],
        "top_name": top["name"],
        "top_R2": top["R2"],
        "top_is_individual": top["name"] == "Individual",
    }


def build_lme_trajectory(rows: list[dict[str, Any]], base_groups: list[str], taxa: list[str]) -> dict[str, Any]:
    """Mean ± 95% CI trajectory per base group × timepoint, per alpha metric.

    Individual patient lines drawn underneath for context.
    Returns a dict keyed by metric name → list of Plotly traces.
    """
    timepoints = sorted_timepoints(rows)
    if len(timepoints) < 2:
        return {}

    result: dict[str, list[dict[str, Any]]] = {}

    for metric, metric_lbl in METRIC_LABELS.items():
        traces: list[dict[str, Any]] = []

        # Individual patient lines (background context)
        patients = sorted(set(r["patient"] for r in rows))
        for bg in base_groups:
            c = base_group_color(bg, base_groups)
            for p in patients:
                p_rows = sorted(
                    [r for r in rows if r["patient"] == p and r.get("base_group", r["group"]) == bg],
                    key=lambda r: r.get("time") or 0,
                )
                if len(p_rows) < 2:
                    continue
                traces.append({
                    "type": "scatter", "mode": "lines+markers", "name": p,
                    "x": [r["timepoint"] for r in p_rows],
                    "y": [metric_value(r, metric, taxa) for r in p_rows],
                    "line": {"color": hex_rgba(c, 0.25), "width": 1},
                    "marker": {"size": 5, "color": hex_rgba(c, 0.25)},
                    "showlegend": False, "hoverinfo": "skip",
                })

        # Group mean ± 95% CI (foreground)
        for bg in base_groups:
            c   = base_group_color(bg, base_groups)
            xs_used, means, uppers, lowers = [], [], [], []
            for tp in timepoints:
                vals = [metric_value(r, metric, taxa) for r in rows
                        if r.get("base_group", r["group"]) == bg and r["timepoint"] == tp]
                if not vals:
                    continue
                m   = float(np.mean(vals))
                sem = float(np.std(vals, ddof=1)) / math.sqrt(len(vals)) if len(vals) > 1 else 0.0
                ci  = 1.96 * sem
                xs_used.append(tp)
                means.append(round(m, 4))
                uppers.append(round(m + ci, 4))
                lowers.append(round(max(0.0, m - ci), 4))

            if not means:
                continue

            # 95% CI band
            traces.append({
                "type": "scatter", "mode": "lines",
                "x": xs_used + xs_used[::-1],
                "y": uppers + lowers[::-1],
                "fill": "toself", "fillcolor": hex_rgba(c, 0.12),
                "line": {"width": 0}, "showlegend": False, "hoverinfo": "skip", "name": bg,
            })
            # Mean line
            traces.append({
                "type": "scatter", "mode": "lines+markers", "name": bg,
                "x": xs_used, "y": means,
                "line": {"color": c, "width": 3},
                "marker": {"size": 10, "color": c, "symbol": "diamond",
                           "line": {"color": "white", "width": 1.5}},
                "error_y": {
                    "type": "data", "symmetric": False,
                    "array":       [u - m for u, m in zip(uppers, means)],
                    "arrayminus":  [m - lo for m, lo in zip(means, lowers)],
                    "color": hex_rgba(c, 0.5), "thickness": 1.5, "width": 6,
                },
                "hovertemplate": (f"<b>{bg}</b><br>%{{x}}<br>"
                                  f"{metric_lbl}: %{{y:.3f}} ± 95%CI<extra></extra>"),
            })

            # Wilcoxon T0 vs T84 annotation
            if len(timepoints) == 2:
                patients_bg = sorted(set(r["patient"] for r in rows
                                         if r.get("base_group", r["group"]) == bg))
                a_vals, b_vals = [], []
                for p in patients_bg:
                    t0v  = [metric_value(r, metric, taxa) for r in rows
                            if r["patient"] == p and (r.get("time") or 0) == 0]
                    t84v = [metric_value(r, metric, taxa) for r in rows
                            if r["patient"] == p and (r.get("time") or 0) > 0]
                    if t0v and t84v:
                        a_vals.append(t0v[0])
                        b_vals.append(t84v[0])
                p_val = wilcoxon_p(a_vals, b_vals)
                lbl   = sig_label(p_val)
                traces.append({
                    "type": "scatter", "mode": "text", "showlegend": False,
                    "x": [timepoints[len(timepoints) // 2]],
                    "y": [max(uppers) * 1.08 if uppers else 1],
                    "text": [f"{bg}: {lbl}"],
                    "textfont": {"size": 9, "color": c},
                    "hoverinfo": "skip",
                })

        result[metric] = traces

    return result
