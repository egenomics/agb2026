"""charts/orchestrator.py — orchestrates all chart-data computation.

Extracted from renderer.py so the renderer is concerned only with HTML
template rendering.  This module owns the compute_chart_data() function and
the ReportConfig dataclass that controls which sections are built.

The public API is re-exported via charts/__init__.py for backward compatibility.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .alpha import build_all_alpha_metrics, build_multimet_alpha, build_rarefaction
from .beta import build_delta_heatmap, build_dendrogram, build_nmds_plot, build_pcoa_chart
from .clinical import (
    build_clinical_correlation,
    build_clinical_slope,
    build_taxa_clinical_heatmap,
)
from .comparative import (
    build_ancom_style,
    build_corr_matrix,
    build_diff_abundance,
    build_heatmap,
    build_volcano,
)
from .distances import bray_curtis_matrix, rows_to_ab
from .individual import (
    build_diversity_rank,
    build_faceted_composition,
    build_nmds_trajectories,
    build_paired_slope,
    build_patient_radar_profiles,
    build_stability_bar,
)
from .insights import generate_chart_explanations, generate_dynamic_insights
from .preprocessing import get_base_groups
from .stats import build_lme_trajectory, build_longitudinal, build_permanova_table, build_stats_table
from .taxonomy import build_donut, build_sunburst, build_taxonomy_views, build_top_taxa


# ── Configuration dataclass ───────────────────────────────────────────────────

_ALL_SECTIONS = frozenset({
    "taxonomy", "alpha", "beta", "individual", "comparative", "stats",
})

@dataclass
class ReportConfig:
    """Controls which chart sections are computed.

    Pass an instance to compute_chart_data() to skip sections you don't need.
    Omit it (or pass None) for the full report — backward-compatible default.

    Example — taxonomy + alpha only::

        cfg = ReportConfig(sections=["taxonomy", "alpha"])
        chart_data = compute_chart_data(result, config=cfg)

    The ``clinical`` section is always auto-enabled when the dataset contains
    clinical fields (sixmwt / il18) regardless of this setting.

    Attributes:
        sections: Which chart sections to compute.  Defaults to all sections.
        max_taxa: Truncate taxa list to this many before chart building.
                  None = use all taxa.  Useful for very large feature tables.
    """
    sections: list[str] = field(default_factory=lambda: sorted(_ALL_SECTIONS))
    max_taxa: int | None = None

    def includes(self, section: str) -> bool:
        return section in self.sections


# ── Main orchestrator ─────────────────────────────────────────────────────────

def compute_chart_data(
    result: Any,
    config: ReportConfig | None = None,
) -> dict[str, Any]:
    """Compute all chart payloads from an IntegrateResult.

    Args:
        result: An IntegrateResult produced by parsers.integrate().
        config: Optional ReportConfig controlling which sections to build.
                Defaults to the full report.

    Returns:
        A dict of chart payloads consumed by render_html().
    """
    cfg         = config or ReportConfig()
    rows        = [r.model_dump() for r in result.rows]
    taxa        = result.taxa[:cfg.max_taxa] if cfg.max_taxa else result.taxa
    groups      = result.groups
    base_groups = get_base_groups(rows)

    data: dict[str, Any] = {
        "meta": {
            "n_samples":    result.n_samples,
            "n_taxa":       result.n_taxa,
            "groups":       groups,
            "base_groups":  base_groups,
            "has_clinical": result.has_clinical,
            "warnings":     result.warnings,
        },
    }

    # Pre-compute Bray-Curtis matrix once — shared by beta, individual, and stats sections
    _need_bc = taxa and (cfg.includes("beta") or cfg.includes("stats") or cfg.includes("individual"))
    _bc_mat = bray_curtis_matrix(rows_to_ab(rows, taxa)) if _need_bc else None

    # ── Taxonomy ──────────────────────────────────────────────────────────────
    if cfg.includes("taxonomy"):
        data["taxonomy_views"] = build_taxonomy_views(rows, taxa, base_groups)
        data["top_taxa"]       = build_top_taxa(rows, taxa)
        data["donut"]          = build_donut(rows, taxa, groups)
        data["sunburst"]       = build_sunburst(rows, taxa, groups)

    # ── Alpha diversity ───────────────────────────────────────────────────────
    if cfg.includes("alpha"):
        data["alpha_metrics"]  = build_all_alpha_metrics(rows, groups, taxa, base_groups)
        data["rarefaction"]    = build_rarefaction(rows, taxa, groups)
        data["multimet_alpha"] = build_multimet_alpha(rows, taxa, groups)

    # ── Beta diversity ────────────────────────────────────────────────────────
    if cfg.includes("beta"):
        bray_traces, bray_pct1, bray_pct2 = build_pcoa_chart(rows, taxa, groups, "bray", _bc_mat)
        jacc_traces, jacc_pct1, jacc_pct2 = build_pcoa_chart(rows, taxa, groups, "jaccard")
        nmds_sa_traces, nmds_sa_pct1, nmds_sa_pct2 = build_nmds_plot(rows, taxa, groups, _bc_mat)

        data["pcoa_bray"]     = {"traces": bray_traces, "pct1": bray_pct1, "pct2": bray_pct2}
        data["pcoa_jaccard"]  = {"traces": jacc_traces, "pct1": jacc_pct1, "pct2": jacc_pct2}
        data["nmds"]          = {"traces": nmds_sa_traces, "pct1": nmds_sa_pct1, "pct2": nmds_sa_pct2}
        data["dendrogram"]    = build_dendrogram(rows, taxa, groups, _bc_mat)
        data["delta_heatmap"] = build_delta_heatmap(rows, taxa)

    # ── Individual / per-patient ──────────────────────────────────────────────
    if cfg.includes("individual"):
        data["paired_slope"]           = build_paired_slope(rows, groups)
        data["stability_bar"]          = build_stability_bar(rows, taxa)
        data["diversity_rank"]         = build_diversity_rank(rows)
        data["patient_radar_profiles"] = build_patient_radar_profiles(rows, taxa)
        data["faceted_composition"]    = build_faceted_composition(rows, taxa)
        _t, _p1, _p2 = build_nmds_trajectories(rows, taxa, _bc_mat)
        data["nmds_trajectories"] = {"traces": _t, "pct1": _p1, "pct2": _p2}

    # ── Comparative ───────────────────────────────────────────────────────────
    if cfg.includes("comparative"):
        data["diff_abundance"] = build_diff_abundance(rows, taxa)
        data["volcano"]        = build_volcano(rows, taxa)
        data["ancom_style"]    = build_ancom_style(rows, taxa)
        data["heatmap"]        = build_heatmap(rows, taxa)
        data["corr_matrix"]    = build_corr_matrix(rows, taxa)

    # ── Longitudinal + statistics ─────────────────────────────────────────────
    if cfg.includes("stats"):
        data["longitudinal"]   = build_longitudinal(rows)
        data["lme_trajectory"] = build_lme_trajectory(rows, base_groups, taxa)
        data["stats_table"]    = build_stats_table(rows, groups, taxa)
        data["permanova"]      = build_permanova_table(rows, taxa, _bc_mat)

    # ── Clinical (auto-enabled when data has clinical fields) ─────────────────
    if result.has_clinical:
        corr_mwt_traces,  r_mwt,  p_mwt  = build_clinical_correlation(rows, "sixmwt")
        corr_il18_traces, r_il18, p_il18  = build_clinical_correlation(rows, "il18")
        data["clinical_sixmwt"]       = build_clinical_slope(rows, "sixmwt")
        data["clinical_il18"]         = build_clinical_slope(rows, "il18")
        data["corr_mwt"]              = {"traces": corr_mwt_traces,  "r": r_mwt,  "p": p_mwt}
        data["corr_il18"]             = {"traces": corr_il18_traces, "r": r_il18, "p": p_il18}
        data["taxa_clinical_heatmap"] = build_taxa_clinical_heatmap(rows, taxa)

    # ── Dynamic insights + per-chart explanations (must run last) ────────────
    data["insights"]      = generate_dynamic_insights(result, data, rows)
    data["explanations"]  = generate_chart_explanations(result, data, rows)

    return data
