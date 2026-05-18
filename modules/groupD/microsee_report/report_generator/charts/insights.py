"""charts/insights.py — dynamic textual insights derived from chart payloads.

Extracted from renderer.py so the renderer is a pure template-filler with no
analytical logic.  All numbers in the returned strings come from real data —
no hard-coded values.
"""
from __future__ import annotations
from typing import Any

import numpy as np

from .preprocessing import get_patient_timepoints, get_base_groups, get_unique_patients, sorted_timepoints


def generate_dynamic_insights(
    result: Any,
    chart_data: dict[str, Any],
    rows: list[dict],
) -> dict[str, str]:
    """Read computed chart payloads and return one human-readable insight per section."""
    taxa = result.taxa
    insights: dict[str, str] = {}

    # ── Taxonomy ──────────────────────────────────────────────────────────────
    if taxa:
        means = {t: float(np.mean([float(r.get(t) or 0) for r in rows])) for t in taxa}
        top_t = max(means, key=lambda t: means[t])
        insights["taxonomy"] = (
            f"{result.n_taxa} families detected across {result.n_samples} samples. "
            f"The most abundant is {top_t} ({means[top_t]:.1f}% mean relative abundance)."
        )

    # ── Alpha diversity ───────────────────────────────────────────────────────
    base_groups = get_base_groups(rows)
    if len(base_groups) >= 2:
        bg_means: dict[str, float] = {}
        for bg in base_groups:
            vals = [float(r.get("shannon") or 0)
                    for r in rows if r.get("base_group", r["group"]) == bg]
            if vals:
                bg_means[bg] = float(np.mean(vals))
        if len(bg_means) >= 2:
            sorted_bgs = sorted(bg_means, key=lambda bg: bg_means[bg], reverse=True)
            a, b = sorted_bgs[0], sorted_bgs[1]
            diff = bg_means[a] - bg_means[b]
            insights["alpha"] = (
                f"Mean Shannon H′: {a}={bg_means[a]:.2f} vs {b}={bg_means[b]:.2f} "
                f"(difference={diff:+.2f}). {a} shows higher alpha diversity."
            )

    # ── Beta diversity / PCoA ─────────────────────────────────────────────────
    bray = chart_data.get("pcoa_bray", {})
    pm   = chart_data.get("permanova", {})
    if bray.get("pct1") is not None:
        cluster_driver = ("individual identity" if pm.get("top_is_individual")
                          else pm.get("top_name", "group membership"))
        insights["pcoa"] = (
            f"PC1 explains {bray['pct1']:.1f}% of variance, "
            f"PC2 {bray.get('pct2', 0):.1f}%. "
            f"Samples cluster primarily by {cluster_driver}."
        )

    # ── PERMANOVA ─────────────────────────────────────────────────────────────
    if pm.get("rows"):
        top_name = pm.get("top_name", "Unknown")
        top_R2   = pm.get("top_R2", 0.0)
        top_row  = next((r for r in pm["rows"] if r[0] == top_name), None)
        if top_row:
            top_p    = float(top_row[3])
            supp_row = next((r for r in pm["rows"] if "Group" in r[0] or "Suppl" in r[0]), None)
            supp_p   = float(supp_row[3]) if supp_row else 1.0
            sig_text = ("Supplementation had no significant effect on community composition."
                        if supp_p > 0.05 else
                        "Supplementation significantly shifted community composition.")
            insights["permanova"] = (
                f"{top_name} explains the most variance "
                f"(R²={top_R2:.3f}, p={top_p:.3f}). {sig_text}"
            )

    # ── Stability ─────────────────────────────────────────────────────────────
    stab = chart_data.get("stability_bar", [])
    if stab and stab[0].get("x") and stab[0].get("y"):
        bc_vals  = stab[0]["x"]
        pts      = stab[0]["y"]
        n_stable = sum(1 for v in bc_vals if float(v) < 0.2)
        median_bc = float(sorted(bc_vals)[len(bc_vals) // 2])
        insights["stability"] = (
            f"Most stable patient: {pts[0]} (BC={float(bc_vals[0]):.3f}). "
            f"Least stable: {pts[-1]} (BC={float(bc_vals[-1]):.3f}). "
            f"Median dissimilarity: {median_bc:.3f}. "
            f"{n_stable}/{len(bc_vals)} patients showed stable composition (BC < 0.2)."
        )

    # ── Delta heatmap ─────────────────────────────────────────────────────────
    patients = get_unique_patients(rows)
    if taxa and patients:
        max_delta   = 0.0
        max_taxon   = ""
        max_patient = ""
        n_stable_delta = 0
        for p in patients:
            r0, r84 = get_patient_timepoints(rows, p)
            if r0 is None or r84 is None:
                continue
            tot0  = sum(float(r0.get(t) or 0) for t in taxa) or 1.0
            tot84 = sum(float(r84.get(t) or 0) for t in taxa) or 1.0
            deltas = {
                t: float(r84.get(t) or 0) / tot84 * 100
                   - float(r0.get(t) or 0) / tot0 * 100
                for t in taxa
            }
            if max(abs(v) for v in deltas.values()) < 5:
                n_stable_delta += 1
            for t, d in deltas.items():
                if abs(d) > abs(max_delta):
                    max_delta   = d
                    max_taxon   = t
                    max_patient = p
        if max_taxon:
            insights["delta"] = (
                f"Most changed family: {max_taxon} (Δ={max_delta:+.1f}% in {max_patient}). "
                f"{n_stable_delta}/{len(patients)} patients showed stable composition "
                f"(all |Δ| < 5%)."
            )

    # ── Longitudinal ─────────────────────────────────────────────────────────
    timepoints = sorted_timepoints(rows)
    if len(timepoints) == 2 and base_groups:
        parts: list[str] = []
        for bg in base_groups:
            t0_vals  = [float(r.get("shannon") or 0) for r in rows
                        if r.get("base_group", r["group"]) == bg
                        and r["timepoint"] == timepoints[0]]
            t84_vals = [float(r.get("shannon") or 0) for r in rows
                        if r.get("base_group", r["group"]) == bg
                        and r["timepoint"] == timepoints[1]]
            if t0_vals and t84_vals:
                delta     = float(np.mean(t84_vals)) - float(np.mean(t0_vals))
                direction = ("increased" if delta > 0.05
                             else "decreased" if delta < -0.05
                             else "remained stable")
                parts.append(f"{bg} {direction} by {abs(delta):.2f}")
        if parts:
            insights["longitudinal"] = (
                "Shannon diversity over time: " + "; ".join(parts) + "."
            )

    # ── Clinical ──────────────────────────────────────────────────────────────
    if result.has_clinical:
        corr_mwt  = chart_data.get("corr_mwt",  {})
        corr_il18 = chart_data.get("corr_il18", {})
        parts_c: list[str] = []
        if corr_mwt.get("r") is not None:
            p_mwt = corr_mwt.get("p") or 1.0
            sig   = "significantly" if float(p_mwt) < 0.05 else "not significantly"
            parts_c.append(
                f"Shannon H′ {sig} correlates with 6MWT "
                f"(r={corr_mwt['r']:.2f}, p={float(p_mwt):.3f})"
            )
        if corr_il18.get("r") is not None:
            p_il18 = corr_il18.get("p") or 1.0
            sig    = "significantly" if float(p_il18) < 0.05 else "not significantly"
            parts_c.append(
                f"and IL-18 (r={corr_il18['r']:.2f}, p={float(p_il18):.3f})"
            )
        if parts_c:
            insights["clinical"] = " ".join(parts_c) + "."

    return insights
