"""charts/insights.py — dynamic textual insights for the HTML report.

Two public functions:
  generate_dynamic_insights   — one short sentence per section (shown in .sec-insight banners)
  generate_chart_explanations — detailed what/finding/pills per chart (shown in ℹ panels)

generate_chart_explanations lives in insights_charts.py and is re-exported here
so callers only need to import from this module.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .insights_charts import generate_chart_explanations
from .preprocessing import (
    get_base_groups,
    get_patient_timepoints,
    get_unique_patients,
    sorted_timepoints,
)

__all__ = ["generate_chart_explanations", "generate_dynamic_insights"]


# ── Section-level one-liners ──────────────────────────────────────────────────


def generate_dynamic_insights(
    result: Any,
    chart_data: dict[str, Any],
    rows: list[dict[str, Any]],
) -> dict[str, str]:
    """Return one human-readable insight per section (used in .sec-insight banners)."""
    taxa = result.taxa
    insights: dict[str, str] = {}

    # Taxonomy
    if taxa:
        means = {t: float(np.mean([float(r.get(t) or 0) for r in rows])) for t in taxa}
        top_t = max(means, key=lambda k: means[k])
        insights["taxonomy"] = (
            f"{result.n_taxa} families detected across {result.n_samples} samples. "
            f"The most abundant is {top_t} ({means[top_t]:.1f}% mean relative abundance)."
        )

    # Alpha diversity
    base_groups = get_base_groups(rows)
    if len(base_groups) >= 2:
        bg_means: dict[str, float] = {}
        for bg in base_groups:
            vals = [
                float(r.get("shannon") or 0)
                for r in rows
                if r.get("base_group", r.get("group")) == bg
            ]
            if vals:
                bg_means[bg] = float(np.mean(vals))
        if len(bg_means) >= 2:
            sorted_bgs = sorted(bg_means, key=lambda k: bg_means[k], reverse=True)
            hi, lo = sorted_bgs[0], sorted_bgs[1]
            diff = bg_means[hi] - bg_means[lo]
            insights["alpha"] = (
                f"Mean Shannon H′: {hi}={bg_means[hi]:.2f} vs {lo}={bg_means[lo]:.2f} "
                f"(difference={diff:+.2f}). {hi} shows higher alpha diversity."
            )

    # Beta diversity / PCoA
    # Use `or {}` (not `.get(key, {})`) so the annotated type stays dict[str, Any].
    bray: dict[str, Any] = chart_data.get("pcoa_bray") or {}
    pm: dict[str, Any] = chart_data.get("permanova") or {}
    if bray.get("pct1") is not None:
        pct1: float = float(bray["pct1"])
        pct2: float = float(bray.get("pct2") or 0)
        cluster_driver: str = (
            "individual identity"
            if pm.get("top_is_individual")
            else str(pm.get("top_name") or "group membership")
        )
        insights["pcoa"] = (
            f"PC1 explains {pct1:.1f}% of variance, "
            f"PC2 {pct2:.1f}%. "
            f"Samples cluster primarily by {cluster_driver}."
        )

    # PERMANOVA
    pm_rows: list[Any] = pm.get("rows") or []
    if pm_rows:
        top_name: str = str(pm.get("top_name") or "Unknown")
        top_R2: float = float(pm.get("top_R2") or 0.0)
        top_row = next((r for r in pm_rows if str(r[0]) == top_name), None)
        if top_row is not None:
            top_p: float = float(top_row[3])
            supp_row = next(
                (r for r in pm_rows if "Group" in str(r[0]) or "Suppl" in str(r[0])),
                None,
            )
            supp_p: float = float(supp_row[3]) if supp_row is not None else 1.0
            sig_text = (
                "Supplementation had no significant effect on community composition."
                if supp_p > 0.05
                else "Supplementation significantly shifted community composition."
            )
            insights["permanova"] = (
                f"{top_name} explains the most variance "
                f"(R²={top_R2:.3f}, p={top_p:.3f}). {sig_text}"
            )

    # Stability
    stab: list[Any] = chart_data.get("stability_bar") or []
    if stab and stab[0].get("x") and stab[0].get("y"):
        bc_vals: list[Any] = stab[0]["x"]
        pts: list[Any] = stab[0]["y"]
        n_stable = sum(1 for v in bc_vals if float(v) < 0.2)
        median_bc = float(sorted(float(v) for v in bc_vals)[len(bc_vals) // 2])
        insights["stability"] = (
            f"Most stable patient: {pts[0]} (BC={float(bc_vals[0]):.3f}). "
            f"Least stable: {pts[-1]} (BC={float(bc_vals[-1]):.3f}). "
            f"Median dissimilarity: {median_bc:.3f}. "
            f"{n_stable}/{len(bc_vals)} patients showed stable composition (BC < 0.2)."
        )

    # Delta heatmap
    patients = get_unique_patients(rows)
    if taxa and patients:
        max_delta, max_taxon, max_patient = 0.0, "", ""
        n_stable_delta = 0
        for p in patients:
            r0, r84 = get_patient_timepoints(rows, p)
            if r0 is None or r84 is None:
                continue
            tot0 = sum(float(r0.get(t) or 0) for t in taxa) or 1.0
            tot84 = sum(float(r84.get(t) or 0) for t in taxa) or 1.0
            deltas = {
                t: float(r84.get(t) or 0) / tot84 * 100 - float(r0.get(t) or 0) / tot0 * 100
                for t in taxa
            }
            if max(abs(v) for v in deltas.values()) < 5:
                n_stable_delta += 1
            for t, d in deltas.items():
                if abs(d) > abs(max_delta):
                    max_delta, max_taxon, max_patient = d, t, p
        if max_taxon:
            insights["delta"] = (
                f"Most changed family: {max_taxon} (Δ={max_delta:+.1f}% in {max_patient}). "
                f"{n_stable_delta}/{len(patients)} patients showed stable composition (all |Δ| < 5%)."
            )

    # Longitudinal
    timepoints = sorted_timepoints(rows)
    if len(timepoints) == 2 and base_groups:
        parts: list[str] = []
        for bg in base_groups:
            t0_vals = [
                float(r.get("shannon") or 0)
                for r in rows
                if r.get("base_group", r.get("group")) == bg and r.get("timepoint") == timepoints[0]
            ]
            t84_vals = [
                float(r.get("shannon") or 0)
                for r in rows
                if r.get("base_group", r.get("group")) == bg and r.get("timepoint") == timepoints[1]
            ]
            if t0_vals and t84_vals:
                delta = float(np.mean(t84_vals)) - float(np.mean(t0_vals))
                direction = (
                    "increased"
                    if delta > 0.05
                    else "decreased"
                    if delta < -0.05
                    else "remained stable"
                )
                parts.append(f"{bg} {direction} by {abs(delta):.2f}")
        if parts:
            insights["longitudinal"] = "Shannon diversity over time: " + "; ".join(parts) + "."

    # Clinical
    if result.has_clinical:
        corr_mwt: dict[str, Any] = chart_data.get("corr_mwt") or {}
        corr_il18: dict[str, Any] = chart_data.get("corr_il18") or {}
        parts_c: list[str] = []
        if corr_mwt.get("r") is not None:
            r_mwt: float = float(corr_mwt["r"])
            p_mwt: float = float(corr_mwt.get("p") or 1.0)
            sig_mwt = "significantly" if p_mwt < 0.05 else "not significantly"
            parts_c.append(
                f"Shannon H′ {sig_mwt} correlates with 6MWT (r={r_mwt:.2f}, p={p_mwt:.3f})"
            )
        if corr_il18.get("r") is not None:
            r_il18: float = float(corr_il18["r"])
            p_il18: float = float(corr_il18.get("p") or 1.0)
            sig_il18 = "significantly" if p_il18 < 0.05 else "not significantly"
            parts_c.append(f"and IL-18 {sig_il18} (r={r_il18:.2f}, p={p_il18:.3f})")
        if parts_c:
            insights["clinical"] = " ".join(parts_c) + "."

    return insights
