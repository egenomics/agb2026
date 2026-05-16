"""charts/stats.py — statistical summary builders (longitudinal, diversity table, PERMANOVA,
                       Wilcoxon / Mann-Whitney tests, LME-style trajectory)."""

from __future__ import annotations
import math
import numpy as np
from .utils import _base_group_color, _sorted_timepoints, group_color, hex_rgba
from .distances import rows_to_ab, bray_curtis_matrix


# ── Non-parametric tests (pure Python, no scipy) ─────────────────────────────

def _wilcoxon_p(a: list[float], b: list[float]) -> float:
    """Two-tailed Wilcoxon signed-rank p-value (exact enumeration for n≤20)."""
    diffs = [float(x) - float(y) for x, y in zip(a, b) if float(x) != float(y)]
    n = len(diffs)
    if n < 2:
        return 1.0
    abs_d = [abs(d) for d in diffs]
    sorted_idx = sorted(range(n), key=lambda i: abs_d[i])
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j < n - 1 and abs_d[sorted_idx[j]] == abs_d[sorted_idx[j + 1]]:
            j += 1
        avg = (i + j + 2) / 2.0
        for k in range(i, j + 1):
            ranks[sorted_idx[k]] = avg
        i = j + 1
    w_plus    = sum(r for r, d in zip(ranks, diffs) if d > 0)
    total_w   = n * (n + 1) / 2.0
    w_obs     = min(w_plus, total_w - w_plus)

    if n <= 20:
        cnt = sum(
            1 for mask in range(1 << n)
            if (s := sum(ranks[i] for i in range(n) if (mask >> i) & 1)) <= w_obs
            or s >= total_w - w_obs
        )
        return round(cnt / (1 << n), 4)
    # Normal approximation with continuity correction
    mu    = total_w / 2
    sigma = math.sqrt(n * (n + 1) * (2 * n + 1) / 24)
    z     = (w_obs - mu + 0.5) / sigma
    az    = abs(z)
    t     = 1 / (1 + 0.2316419 * az)
    pd    = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))))
    return round(min(1.0, 2 * pd * math.exp(-0.5 * az ** 2) / math.sqrt(2 * math.pi)), 4)


def _mannwhitney_p(a: list[float], b: list[float]) -> float:
    """Two-tailed Mann-Whitney U p-value (normal approximation)."""
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return 1.0
    combined = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    n        = len(combined)
    ranks    = [0.0] * n
    i        = 0
    while i < n:
        j = i
        while j < n - 1 and combined[j][0] == combined[j + 1][0]:
            j += 1
        avg = (i + j + 2) / 2.0
        for k in range(i, j + 1):
            ranks[k] = avg
        i = j + 1
    r_a   = sum(r for r, (_, g) in zip(ranks, combined) if g == 0)
    u_a   = r_a - na * (na + 1) / 2
    u     = min(u_a, na * nb - u_a)
    mu    = na * nb / 2
    sigma = math.sqrt(na * nb * (na + nb + 1) / 12)
    if sigma == 0:
        return 1.0
    z  = (u - mu) / sigma
    az = abs(z)
    t  = 1 / (1 + 0.2316419 * az)
    pd = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))))
    return round(min(1.0, 2 * pd * math.exp(-0.5 * az ** 2) / math.sqrt(2 * math.pi)), 4)


def _sig_label(p: float) -> str:
    if p < 0.001: return "p<0.001 ***"
    if p < 0.01:  return f"p={p:.3f} **"
    if p < 0.05:  return f"p={p:.3f} *"
    return f"p={p:.3f} ns"


# ── Existing chart builders ───────────────────────────────────────────────────

def build_longitudinal(rows: list[dict]) -> list[dict]:
    """Group mean Shannon per timepoint, connected across time."""
    timepoints  = _sorted_timepoints(rows)
    base_groups = sorted(set(r.get("base_group", r["group"]) for r in rows))
    traces = []
    for bg in base_groups:
        c = _base_group_color(bg, base_groups)
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


def build_stats_table(rows: list[dict], groups: list[str], taxa: list[str] | None = None) -> dict:
    from .alpha import METRIC_LABELS, _metric_val
    metrics = list(METRIC_LABELS.keys())
    header = ["Group", "N"] + [f"{METRIC_LABELS[m]} mean ± SD" for m in metrics]
    table_rows = []
    for g in groups:
        g_rows = [r for r in rows if r["group"] == g]
        cells: list = [g, len(g_rows)]
        for m in metrics:
            vals = [_metric_val(r, m, taxa) for r in g_rows]
            non_zero = [v for v in vals if v > 0]
            if non_zero:
                cells.append(f"{np.mean(non_zero):.3f} ± {np.std(non_zero):.3f}")
            else:
                cells.append("—")
        table_rows.append(cells)
    return {"header": header, "rows": table_rows}


def build_permanova_table(rows: list[dict], taxa: list[str]) -> dict:
    """PERMANOVA: pseudo-F, R², p-value (99 permutations, Bray-Curtis)."""
    n = len(rows)
    if n < 4 or not taxa:
        return {"header": [], "rows": [], "top_name": "", "top_R2": 0.0, "top_is_individual": False}

    ab   = rows_to_ab(rows, taxa)
    mat  = bray_curtis_matrix(ab)
    flat = mat.flatten()

    def pseudo_f(labels: list[str]) -> tuple[float, float]:
        gs    = list(set(labels))
        ss_tot = sum(flat[i * n + j] ** 2 for i in range(n) for j in range(i + 1, n)) / n
        ss_w  = 0.0
        for g in gs:
            idx = [i for i, l in enumerate(labels) if l == g]
            ng  = len(idx)
            if ng < 2:
                continue
            ss_w += sum(flat[idx[a] * n + idx[b]] ** 2
                        for a in range(len(idx)) for b in range(a + 1, len(idx))) / ng
        ss_b = ss_tot - ss_w
        k    = len(gs)
        R2   = ss_b / (ss_tot or 1.0)
        F    = (ss_b / max(k - 1, 1)) / max(ss_w / max(n - k, 1), 1e-10)
        return F, R2

    _rng = [sum(hash(r.get("sample_id", "")) & 0x7fffffff for r in rows) % 0x7fffffff or 42]
    def _rand() -> float:
        _rng[0] = (_rng[0] * 1664525 + 1013904223) & 0x7fffffff
        return _rng[0] / 0x7fffffff

    def permute(labels: list[str], n_perm: int = 99) -> tuple[float, float, float]:
        obs_F, obs_R2 = pseudo_f(labels)
        lab  = list(labels)
        cnt  = 1
        for _ in range(n_perm):
            for j in range(len(lab) - 1, 0, -1):
                k = int(_rand() * (j + 1))
                lab[j], lab[k] = lab[k], lab[j]
            if pseudo_f(lab)[0] >= obs_F:
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


# ── LME-style trajectory ──────────────────────────────────────────────────────

def build_lme_trajectory(rows: list[dict], base_groups: list[str], taxa: list[str]) -> dict:
    """
    Mean ± 95% CI trajectory per base group × timepoint, per alpha metric.
    Individual patient lines drawn underneath (for context).
    Returns a dict keyed by metric name → list of Plotly traces.
    """
    from .alpha import METRIC_LABELS, _metric_val  # local import avoids top-level circular dep

    timepoints = _sorted_timepoints(rows)
    if len(timepoints) < 2:
        return {}

    result: dict[str, list[dict]] = {}

    for metric, metric_lbl in METRIC_LABELS.items():
        traces: list[dict] = []

        # Individual patient lines (background)
        patients = sorted(set(r["patient"] for r in rows))
        for bg in base_groups:
            c = _base_group_color(bg, base_groups)
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
                    "y": [_metric_val(r, metric, taxa) for r in p_rows],
                    "line": {"color": hex_rgba(c, 0.25), "width": 1},
                    "marker": {"size": 5, "color": hex_rgba(c, 0.25)},
                    "showlegend": False, "hoverinfo": "skip",
                })

        # Group mean ± 95% CI (foreground)
        for bg in base_groups:
            c   = _base_group_color(bg, base_groups)
            xs_used, means, uppers, lowers = [], [], [], []
            for tp in timepoints:
                vals = [_metric_val(r, metric, taxa) for r in rows
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
                    "arrayminus":  [m - l for m, l in zip(means, lowers)],
                    "color": hex_rgba(c, 0.5), "thickness": 1.5, "width": 6,
                },
                "hovertemplate": (f"<b>{bg}</b><br>%{{x}}<br>"
                                  f"{metric_lbl}: %{{y:.3f}} ± 95%CI<extra></extra>"),
            })

            # Test: paired Wilcoxon T0 vs T84
            if len(timepoints) == 2:
                patients_bg = sorted(set(r["patient"] for r in rows
                                         if r.get("base_group", r["group"]) == bg))
                a_vals, b_vals = [], []
                for p in patients_bg:
                    t0v  = [_metric_val(r, metric, taxa) for r in rows
                            if r["patient"] == p and (r.get("time") or 0) == 0]
                    t84v = [_metric_val(r, metric, taxa) for r in rows
                            if r["patient"] == p and (r.get("time") or 0) > 0]
                    if t0v and t84v:
                        a_vals.append(t0v[0]); b_vals.append(t84v[0])
                p_val = _wilcoxon_p(a_vals, b_vals)
                lbl   = _sig_label(p_val)
                # Annotation between the two timepoints
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
