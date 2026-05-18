"""charts/comparative.py — comparative chart builders (LFC, volcano, heatmap, correlation matrix)."""
from __future__ import annotations
import math
import numpy as np
from .utils import hex_rgba, base_group_color
from .distances import rows_to_ab
from .stats_helpers import welch_ttest_p, bh_fdr, wilcoxon_p, sig_label


def build_diff_abundance(rows: list[dict], taxa: list[str]) -> list[dict]:
    """Log2 fold change per taxon: T84 vs T0 within each base group."""
    base_groups = sorted(set(r.get("base_group", r["group"]) for r in rows))
    traces = []
    for bg in base_groups:
        t0_rows  = [r for r in rows if r.get("base_group", r["group"]) == bg and r.get("timepoint") == "T0"]
        t84_rows = [r for r in rows if r.get("base_group", r["group"]) == bg and r.get("timepoint") == "T84"]
        if not t0_rows or not t84_rows:
            continue
        stats = [(t,
                  float(np.mean([float(r.get(t) or 0) for r in t0_rows])),
                  float(np.mean([float(r.get(t) or 0) for r in t84_rows]))) for t in taxa]
        lfc_vals = [round(math.log2((m84 + 0.1) / (m0 + 0.1)), 3) for _, m0, m84 in stats]
        traces.append({
            "type": "bar", "orientation": "h", "name": bg,
            "x": lfc_vals,
            "y": [t for t, *_ in stats],
            "marker": {"color": ["#D84E6A" if v > 0 else "#4A7ED4" for v in lfc_vals]},
            "hovertemplate": "%{customdata}<extra></extra>",
            "customdata": [f"<b>{t}</b><br>T0: {m0:.2f}%  T84: {m84:.2f}%<br>Log2FC: {lfc_v:.3f}"
                           for (t, m0, m84), lfc_v in zip(stats, lfc_vals)],
        })
    return traces


def build_volcano(rows: list[dict], taxa: list[str]) -> list[dict]:
    """Volcano: log2FC vs –log10(p); points coloured by BH-FDR q < 0.1."""
    base_groups = sorted(set(r.get("base_group", r["group"]) for r in rows))
    traces = []
    for bg in base_groups:
        t0_rows  = [r for r in rows if r.get("base_group", r["group"]) == bg and r.get("timepoint") == "T0"]
        t84_rows = [r for r in rows if r.get("base_group", r["group"]) == bg and r.get("timepoint") == "T84"]
        if not t0_rows or not t84_rows:
            continue
        c = base_group_color(bg, base_groups)
        raw = []
        for t in taxa:
            a = [float(r.get(t) or 0) for r in t0_rows]
            b = [float(r.get(t) or 0) for r in t84_rows]
            lfc = round(math.log2((float(np.mean(b)) + 0.1) / (float(np.mean(a)) + 0.1)), 3)
            raw.append((t, lfc, welch_ttest_p(a, b)))
        p_vals  = [p for _, _, p in raw]
        q_vals  = bh_fdr(p_vals)
        vstats  = [(t, lfc, round(-math.log10(max(p, 1e-10)), 3), q)
                   for (t, lfc, p), q in zip(raw, q_vals)]
        sig     = [abs(lfc) > 0.5 and q < 0.1 for _, lfc, _, q in vstats]
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
            "hovertemplate": ("<b>%{text}</b><br>Log2FC: %{x:.3f}<br>"
                              "−log10(p): %{y:.2f}<br>FDR q: %{customdata[0]:.3f}<extra></extra>"),
        })
    return traces


def build_ancom_style(rows: list[dict], taxa: list[str]) -> list[dict]:
    """CLR-transformed paired Wilcoxon (T0→T84) per taxon, per base group.

    Centered log-ratio is the compositionally-correct transform for differential
    abundance — equivalent in spirit to ANCOM/ANCOM-BC without the bias correction
    step.  Results displayed as an effect-size bar (CLR mean diff) coloured by FDR.
    """
    base_groups = sorted(set(r.get("base_group", r["group"]) for r in rows))
    traces: list[dict] = []

    def _clr(r: dict) -> dict:
        """CLR-transform a single sample row (pseudocount 0.01)."""
        vals = {t: float(r.get(t) or 0) + 0.01 for t in taxa}
        gm   = math.exp(sum(math.log(v) for v in vals.values()) / len(vals))
        return {t: math.log(v / gm) for t, v in vals.items()}

    for bg in base_groups:
        c       = base_group_color(bg, base_groups)
        bg_rows = [r for r in rows if r.get("base_group", r["group"]) == bg]
        clr_rows = {r["sample_id"]: _clr(r) for r in bg_rows}

        patients = sorted(set(r["patient"] for r in bg_rows))
        taxon_stats = []
        p_vals_raw  = []
        for t in taxa:
            a, b = [], []
            for p in patients:
                t0r  = [clr_rows[r["sample_id"]][t] for r in bg_rows
                        if r["patient"] == p and (r.get("time") or 0) == 0 and r["sample_id"] in clr_rows]
                t84r = [clr_rows[r["sample_id"]][t] for r in bg_rows
                        if r["patient"] == p and (r.get("time") or 0) > 0  and r["sample_id"] in clr_rows]
                if t0r and t84r:
                    a.append(t0r[0]); b.append(t84r[0])
            diff  = round(float(np.mean(b)) - float(np.mean(a)), 4) if a else 0.0
            p_raw = wilcoxon_p(a, b)
            taxon_stats.append((t, diff, p_raw))
            p_vals_raw.append(p_raw)

        q_vals  = bh_fdr(p_vals_raw)
        results = [(t, diff, p, q) for (t, diff, p), q in zip(taxon_stats, q_vals)]

        def _color(diff: float, q: float) -> str:
            if q < 0.1 and diff > 0:  return "#D84E6A"
            if q < 0.1 and diff < 0:  return "#4A7ED4"
            return hex_rgba(c, 0.45)

        traces.append({
            "type": "bar", "orientation": "h", "name": bg,
            "x": [diff for _, diff, _, _ in results],
            "y": [t for t, *_ in results],
            "marker": {"color": [_color(d, q) for _, d, _, q in results]},
            "customdata": [[p, q] for _, _, p, q in results],
            "hovertemplate": (
                "<b>%{y}</b><br>CLR mean diff: %{x:.3f}<br>"
                "p=%{customdata[0]:.3f}  FDR q=%{customdata[1]:.3f}<extra></extra>"
            ),
        })

    return traces


def build_heatmap(rows: list[dict], taxa: list[str]) -> list[dict]:
    return [{
        "type": "heatmap",
        "x": list(taxa),
        "y": [r["sample_id"] for r in rows],
        "z": [[round(float(r.get(t) or 0), 2) for t in taxa] for r in rows],
        "colorscale": "Blues",
        "hovertemplate": "<b>%{y}</b> · <b>%{x}</b><br>%{z:.2f}%<extra></extra>",
    }]


def build_corr_matrix(rows: list[dict], taxa: list[str]) -> list[dict]:
    ab = rows_to_ab(rows, taxa)
    if ab.shape[0] < 2:
        return []
    corr = np.corrcoef(ab.T)
    return [{
        "type": "heatmap",
        "x": list(taxa), "y": list(taxa),
        "z": [[round(float(corr[i, j]), 3) for j in range(len(taxa))] for i in range(len(taxa))],
        "colorscale": "RdBu", "zmid": 0, "zmin": -1, "zmax": 1,
        "hovertemplate": "<b>%{x}</b> vs <b>%{y}</b><br>r = %{z:.3f}<extra></extra>",
    }]
