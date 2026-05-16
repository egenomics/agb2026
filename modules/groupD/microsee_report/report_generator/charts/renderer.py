"""charts/renderer.py — aggregates all chart data and renders the self-contained HTML report."""

from __future__ import annotations
import json
from pathlib import Path
from typing import Any

# Load Plotly.js from the bundled file so the report works offline on HPC nodes.
_PLOTLY_JS_PATH = Path(__file__).parent / "plotly.min.js"
if _PLOTLY_JS_PATH.exists():
    _PLOTLY_JS_INLINE = f"<script>{_PLOTLY_JS_PATH.read_text(encoding='utf-8')}</script>"
else:
    # Attempt to auto-download (requires outbound internet — not available on most HPC nodes).
    _PLOTLY_URL = "https://cdn.plot.ly/plotly-2.35.2.min.js"
    try:
        import urllib.request
        print(f"[MicroSee] plotly.min.js not found — downloading from {_PLOTLY_URL} ...")
        urllib.request.urlretrieve(_PLOTLY_URL, _PLOTLY_JS_PATH)
        _PLOTLY_JS_INLINE = f"<script>{_PLOTLY_JS_PATH.read_text(encoding='utf-8')}</script>"
        print("[MicroSee] plotly.min.js downloaded and cached.")
    except Exception:
        import warnings
        warnings.warn(
            "\n\ncharts/plotly.min.js not found and download failed.\n"
            "On HPC nodes without internet, commit the file to git first:\n"
            "  curl -fsSL https://cdn.plot.ly/plotly-2.35.2.min.js \\\n"
            "       -o report_generator/charts/plotly.min.js\n"
            "  git add report_generator/charts/plotly.min.js\n"
            "Falling back to CDN — charts will be BLANK on offline nodes.\n",
            RuntimeWarning, stacklevel=2,
        )
        _PLOTLY_JS_INLINE = f'<script src="{_PLOTLY_URL}"></script>'

# Load the HTML template from the separate file
_TEMPLATE = (Path(__file__).parent / "template.html").read_text(encoding="utf-8")

from .config import BASE_LAYOUT, BASE_CONFIG, THEME
from .taxonomy   import (build_taxonomy_stacked, build_taxonomy_views,
                         build_top_taxa, build_donut, build_sunburst)
from .alpha      import (build_all_alpha_metrics, build_rarefaction, build_multimet_alpha)
from .beta       import build_pcoa_chart, build_nmds_plot, build_dendrogram, build_delta_heatmap
from .individual import (build_paired_slope, build_stability_bar, build_diversity_rank,
                         build_patient_radar, build_faceted_composition, build_nmds_trajectories)
from .comparative import build_diff_abundance, build_volcano, build_heatmap, build_corr_matrix, build_ancom_style
from .clinical   import build_clinical_slope, build_clinical_correlation, build_taxa_clinical_heatmap
from .stats      import build_longitudinal, build_stats_table, build_permanova_table, build_lme_trajectory


def compute_chart_data(result: Any) -> dict[str, Any]:
    rows        = [r.model_dump() for r in result.rows]
    taxa        = result.taxa
    groups      = result.groups
    base_groups = sorted(set(r.get("base_group", r["group"]) for r in rows))

    bray_traces, bray_pct1, bray_pct2 = build_pcoa_chart(rows, taxa, groups, "bray")
    jacc_traces, jacc_pct1, jacc_pct2 = build_pcoa_chart(rows, taxa, groups, "jaccard")
    nmds_traces, nmds_pct1, nmds_pct2 = build_nmds_trajectories(rows, taxa)
    nmds_sa_traces, nmds_sa_pct1, nmds_sa_pct2 = build_nmds_plot(rows, taxa, groups)

    data: dict[str, Any] = {
        "meta": {
            "n_samples":    result.n_samples,
            "n_taxa":       result.n_taxa,
            "groups":       groups,
            "base_groups":  base_groups,
            "has_clinical": result.has_clinical,
            "warnings":     result.warnings,
        },
        # Taxonomy — static charts + all filter variants for the stacked bar
        "taxonomy_views":   build_taxonomy_views(rows, taxa, base_groups),
        "top_taxa":         build_top_taxa(rows, taxa),
        "donut":            build_donut(rows, taxa, groups),
        "sunburst":         build_sunburst(rows, taxa, groups),
        # Alpha diversity — all metrics pre-computed for JS toggle (includes Wilcoxon brackets)
        "alpha_metrics":    build_all_alpha_metrics(rows, groups, taxa, base_groups),
        "rarefaction":      build_rarefaction(rows, taxa, groups),
        "multimet_alpha":   build_multimet_alpha(rows, taxa, groups),
        # Beta diversity
        "pcoa_bray":        {"traces": bray_traces, "pct1": bray_pct1, "pct2": bray_pct2},
        "pcoa_jaccard":     {"traces": jacc_traces, "pct1": jacc_pct1, "pct2": jacc_pct2},
        "nmds":             {"traces": nmds_sa_traces, "pct1": nmds_sa_pct1, "pct2": nmds_sa_pct2},
        "dendrogram":       build_dendrogram(rows, taxa, groups),
        "delta_heatmap":    build_delta_heatmap(rows, taxa),
        # Individual
        "paired_slope":        build_paired_slope(rows, groups),
        "stability_bar":       build_stability_bar(rows, taxa),
        "diversity_rank":      build_diversity_rank(rows),
        "patient_radar":       build_patient_radar(rows, taxa),
        "faceted_composition": build_faceted_composition(rows, taxa),
        "nmds_trajectories":   {"traces": nmds_traces, "pct1": nmds_pct1, "pct2": nmds_pct2},
        # Comparative
        "diff_abundance":   build_diff_abundance(rows, taxa),
        "volcano":          build_volcano(rows, taxa),
        "ancom_style":      build_ancom_style(rows, taxa),
        "heatmap":          build_heatmap(rows, taxa),
        "corr_matrix":      build_corr_matrix(rows, taxa),
        # Longitudinal + statistics
        "longitudinal":     build_longitudinal(rows),
        "lme_trajectory":   build_lme_trajectory(rows, base_groups, taxa),
        "stats_table":      build_stats_table(rows, groups, taxa),
        "permanova":        build_permanova_table(rows, taxa),
    }

    if result.has_clinical:
        corr_mwt_traces,  r_mwt,  p_mwt  = build_clinical_correlation(rows, "sixmwt")
        corr_il18_traces, r_il18, p_il18  = build_clinical_correlation(rows, "il18")
        data["clinical_sixmwt"]        = build_clinical_slope(rows, "sixmwt")
        data["clinical_il18"]          = build_clinical_slope(rows, "il18")
        data["corr_mwt"]               = {"traces": corr_mwt_traces,  "r": r_mwt,  "p": p_mwt}
        data["corr_il18"]              = {"traces": corr_il18_traces, "r": r_il18, "p": p_il18}
        data["taxa_clinical_heatmap"]  = build_taxa_clinical_heatmap(rows, taxa)

    return data


def render_html(chart_data: dict) -> str:
    has_clinical = chart_data["meta"]["has_clinical"]
    meta         = chart_data["meta"]
    base_groups  = meta.get("base_groups", [])

    clinical_nav = '<a href="#sec-clinical">Clinical</a>' if has_clinical else ""

    clinical_section = ""
    if has_clinical:
        clinical_section = """\
<section id="sec-clinical">
  <div class="sec-header">Clinical Outcomes</div>
  <div class="grid2">
    <div class="chart-card">
      <div class="chart-title">6-Min Walk Test</div>
      <div class="chart-sub">Per patient change · group mean dashed</div>
      <div id="chart-mwt" class="plot-div"></div>
    </div>
    <div class="chart-card">
      <div class="chart-title">IL-18 Cytokine</div>
      <div class="chart-sub">Per patient change · group mean dashed</div>
      <div id="chart-il18" class="plot-div"></div>
    </div>
  </div>
  <div class="grid2">
    <div class="chart-card" id="card-corr-mwt">
      <div class="chart-title">Shannon H′ vs 6MWT</div>
      <div class="chart-sub" id="sub-corr-mwt">Pearson correlation · dashed = regression line</div>
      <div id="chart-corr-mwt" class="plot-div"></div>
    </div>
    <div class="chart-card" id="card-corr-il18">
      <div class="chart-title">Shannon H′ vs IL-18</div>
      <div class="chart-sub" id="sub-corr-il18">Pearson correlation · dashed = regression line</div>
      <div id="chart-corr-il18" class="plot-div"></div>
    </div>
  </div>
  <div class="grid1">
    <div class="chart-card">
      <div class="chart-title">Taxa × Clinical Correlation Heatmap</div>
      <div class="chart-sub">Spearman ρ between Δ taxon abundance and Δ clinical outcome (T84 − T0) · * p&lt;0.05 · ** p&lt;0.01</div>
      <div id="chart-taxa-clinical" class="plot-div-lg"></div>
    </div>
  </div>
</section>"""

    warnings_html = "".join(f'<div class="warning">{w}</div>' for w in meta["warnings"])

    group_filter_buttons = "\n".join(
        f'<button class="ctrl-btn grp-btn" data-grp="{bg}">{bg}</button>'
        for bg in base_groups
    )

    html = _TEMPLATE
    for placeholder, value in [
        ("__FONT__",                 THEME["font"]),
        ("__BG__",                   THEME["bg"]),
        ("__PAPER__",                THEME["paper"]),
        ("__TEXT__",                 THEME["text"]),
        ("__TEXT2__",                THEME["text2"]),
        ("__N_SAMPLES__",            str(meta["n_samples"])),
        ("__N_TAXA__",               str(meta["n_taxa"])),
        ("__N_GROUPS__",             str(len(meta["groups"]))),
        ("__GROUPS_STR__",           ", ".join(meta["groups"])),
        ("__WARNINGS__",             warnings_html),
        ("__CLINICAL_NAV__",         clinical_nav),
        ("__CLINICAL_SECTION__",     clinical_section),
        ("__GROUP_FILTER_BUTTONS__", group_filter_buttons),
        ("__DATA_JSON__",            json.dumps(chart_data,  allow_nan=False)),
        ("__LAYOUT_JSON__",          json.dumps(BASE_LAYOUT, allow_nan=False)),
        ("__CONFIG_JSON__",          json.dumps(BASE_CONFIG,  allow_nan=False)),
        ("__PLOTLY_SCRIPT__",        _PLOTLY_JS_INLINE),
    ]:
        html = html.replace(placeholder, value)
    return html
