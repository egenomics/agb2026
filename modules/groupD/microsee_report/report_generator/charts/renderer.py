"""charts/renderer.py — aggregates all chart data and renders the self-contained HTML report."""

from __future__ import annotations
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

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


def generate_dynamic_insights(result: Any, chart_data: dict[str, Any]) -> dict[str, str]:
    """Read computed data and return one human-readable insight string per chart section.

    All numbers in the returned strings come from real data — no hard-coded values.
    """
    rows = [r.model_dump() for r in result.rows]
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
    base_groups = sorted(set(r.get("base_group", r["group"]) for r in rows))
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

    # ── Beta diversity / PCoA ────────────────────────────────────────────────
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

    # ── PERMANOVA ────────────────────────────────────────────────────────────
    if pm.get("rows"):
        top_name = pm.get("top_name", "Unknown")
        top_R2   = pm.get("top_R2", 0.0)
        top_row  = next((r for r in pm["rows"] if r[0] == top_name), None)
        if top_row:
            top_p = float(top_row[3])
            supp_row = next((r for r in pm["rows"] if "Group" in r[0] or "Suppl" in r[0]), None)
            supp_p = float(supp_row[3]) if supp_row else 1.0
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
        bc_vals = stab[0]["x"]
        pts     = stab[0]["y"]
        n_stable = sum(1 for v in bc_vals if float(v) < 0.2)
        median_bc = float(sorted(bc_vals)[len(bc_vals) // 2])
        insights["stability"] = (
            f"Most stable patient: {pts[0]} (BC={float(bc_vals[0]):.3f}). "
            f"Least stable: {pts[-1]} (BC={float(bc_vals[-1]):.3f}). "
            f"Median dissimilarity: {median_bc:.3f}. "
            f"{n_stable}/{len(bc_vals)} patients showed stable composition (BC < 0.2)."
        )

    # ── Delta heatmap ─────────────────────────────────────────────────────────
    patients = sorted(set(r["patient"] for r in rows))
    if taxa and patients:
        max_delta = 0.0
        max_taxon = ""
        max_patient = ""
        n_stable_delta = 0
        for p in patients:
            r0s  = [r for r in rows if r["patient"] == p and (r.get("time") or 0) == 0]
            r84s = sorted([r for r in rows if r["patient"] == p and (r.get("time") or 0) > 0],
                          key=lambda r: r.get("time") or 0, reverse=True)
            if not r0s or not r84s:
                continue
            r0, r84 = r0s[0], r84s[0]
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
                    max_delta = d
                    max_taxon = t
                    max_patient = p
        if max_taxon:
            insights["delta"] = (
                f"Most changed family: {max_taxon} (Δ={max_delta:+.1f}% in {max_patient}). "
                f"{n_stable_delta}/{len(patients)} patients showed stable composition "
                f"(all |Δ| < 5%)."
            )

    # ── Longitudinal ─────────────────────────────────────────────────────────
    timepoints = sorted(
        set(r["timepoint"] for r in rows),
        key=lambda tp: min(r.get("time") or 0 for r in rows if r["timepoint"] == tp),
    )
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
                delta = float(np.mean(t84_vals)) - float(np.mean(t0_vals))
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


def _patient_bc(rows: list[dict], taxa: list[str], patient: str) -> float:
    """Compute Bray-Curtis dissimilarity between T0 and T84 for one patient."""
    r0s  = [r for r in rows if r["patient"] == patient and (r.get("time") or 0) == 0]
    r84s = sorted(
        [r for r in rows if r["patient"] == patient and (r.get("time") or 0) > 0],
        key=lambda r: r.get("time") or 0, reverse=True,
    )
    if not r0s or not r84s or not taxa:
        return 0.0
    ab0  = np.array([float(r0s[0].get(t) or 0) for t in taxa])
    ab84 = np.array([float(r84s[0].get(t) or 0) for t in taxa])
    denom = float(np.sum(ab0 + ab84))
    return round(float(np.sum(np.abs(ab0 - ab84))) / denom, 3) if denom > 0 else 0.0


def render_patient_html(patient_id: str, result: Any) -> str:
    """Generate a self-contained per-patient HTML report.

    Contains: stability score, radar T0 vs T84 vs group mean, composition bars,
    clinical outcomes (if available), and a plain-language interpretation.
    """
    from .individual import build_patient_radar, build_faceted_composition
    from .clinical   import build_clinical_slope

    rows  = [r.model_dump() for r in result.rows]
    taxa  = result.taxa
    p_rows = [r for r in rows if r["patient"] == patient_id]
    if not p_rows:
        raise ValueError(f"Patient {patient_id!r} not found in result")

    bc    = _patient_bc(rows, taxa, patient_id)
    group = p_rows[0].get("base_group", p_rows[0]["group"])
    stab_label = (
        "very stable"          if bc < 0.10 else
        "stable"               if bc < 0.20 else
        "moderately variable"  if bc < 0.35 else
        "highly variable"
    )

    # Build charts — radar uses all patients so group mean is correct
    radar_traces = build_patient_radar(rows, taxa)
    comp_data    = build_faceted_composition(p_rows, taxa)

    patient_chart_data: dict[str, Any] = {
        "radar": radar_traces,
        "comp":  comp_data,
        "bc":    bc,
    }
    if result.has_clinical:
        patient_chart_data["sixmwt"] = build_clinical_slope(
            [r for r in rows if r.get("sixmwt") is not None], "sixmwt"
        )
        patient_chart_data["il18"] = build_clinical_slope(
            [r for r in rows if r.get("il18") is not None], "il18"
        )

    layout_json  = json.dumps(BASE_LAYOUT, allow_nan=False)
    config_json  = json.dumps(BASE_CONFIG,  allow_nan=False)
    data_json    = json.dumps(patient_chart_data, allow_nan=False)

    clinical_block = ""
    if result.has_clinical:
        clinical_block = """\
      <div class="grid2">
        <div class="chart-card">
          <div class="chart-title">6-Min Walk Test</div>
          <div id="pt-mwt" class="plot-div"></div>
        </div>
        <div class="chart-card">
          <div class="chart-title">IL-18 Cytokine</div>
          <div id="pt-il18" class="plot-div"></div>
        </div>
      </div>"""

    clinical_js = ""
    if result.has_clinical:
        clinical_js = """\
    if(D.sixmwt&&D.sixmwt.length)
      Plotly.newPlot('pt-mwt', D.sixmwt, merge(L,{height:320,yaxis:{title:{text:'6MWT (m)'}}}), C);
    if(D.il18&&D.il18.length)
      Plotly.newPlot('pt-il18', D.il18, merge(L,{height:320,yaxis:{title:{text:'IL-18 (pg/mL)'}}}), C);"""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>MicroSee — Patient {patient_id}</title>
{_PLOTLY_JS_INLINE}
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:{BASE_LAYOUT['font']['family']};background:#FEF3EC;color:#6B3A2A;padding:24px}}
h1{{font-size:22px;font-weight:800;color:#3E1A0E;margin-bottom:4px}}
.subtitle{{font-size:13px;color:#8B5860;margin-bottom:24px}}
.score-card{{background:#fff;border-radius:12px;padding:20px 24px;margin-bottom:20px;
  box-shadow:0 1px 6px rgba(107,58,42,.10);display:flex;align-items:center;gap:20px}}
.score-num{{font-size:48px;font-weight:800;color:#D97A3A;line-height:1}}
.score-right{{display:flex;flex-direction:column;gap:4px}}
.score-label{{font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:#8B5860}}
.score-text{{font-size:13px;color:#6B3A2A;max-width:460px;line-height:1.5}}
.grid2{{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin-bottom:16px}}
.chart-card{{background:#fff;border-radius:10px;padding:16px;
  box-shadow:0 1px 4px rgba(107,58,42,.08)}}
.chart-title{{font-size:13px;font-weight:700;color:#3E1A0E;margin-bottom:8px}}
.plot-div{{width:100%;height:340px}}
.insight{{background:#FFF8F3;border-left:4px solid #D97A3A;border-radius:0 8px 8px 0;
  padding:12px 16px;font-size:13px;line-height:1.6;margin-top:16px}}
</style>
</head>
<body>
<h1>Patient {patient_id}</h1>
<div class="subtitle">Group: {group}</div>

<div class="score-card">
  <div class="score-num">{bc:.3f}</div>
  <div class="score-right">
    <span class="score-label">Microbiome Stability (Bray-Curtis T0 → T84)</span>
    <span class="score-text">
      Your microbiome was <strong>{stab_label}</strong> between visits
      (BC={bc:.3f} — 0 = identical, 1 = completely different).
      {'A low score indicates consistent microbial community structure.' if bc < 0.2
       else 'This degree of change is common in dietary intervention studies.'}
    </span>
  </div>
</div>

<div class="grid2">
  <div class="chart-card">
    <div class="chart-title">Compositional Profile — Radar</div>
    <div id="pt-radar" class="plot-div"></div>
  </div>
  <div class="chart-card">
    <div class="chart-title">Composition T0 vs T84</div>
    <div id="pt-comp" class="plot-div"></div>
  </div>
</div>

{clinical_block}

<div class="insight">
  <strong>Key finding:</strong>
  Patient {patient_id} ({group}) had a microbiome stability score of {bc:.3f}.
  {'Changes were minimal and within normal day-to-day variation.' if bc < 0.2
   else f'The largest compositional shift was detected, warranting closer monitoring.'}
  {'Clinical markers were available for this patient and are shown above.' if result.has_clinical else ''}
</div>

<script>
(function(){{
  var D={data_json};
  var L={layout_json};
  var C={config_json};
  function merge(b,e){{return Object.assign(JSON.parse(JSON.stringify(b)),e);}}
  if(D.radar&&D.radar.length)
    Plotly.newPlot('pt-radar',
      D.radar.filter(function(t){{return t.name&&t.name.includes('{patient_id}')||t.name&&t.name.includes('{group}');}}),
      merge(L,{{height:340,polar:{{radialaxis:{{visible:true,ticksuffix:'%'}}}},showlegend:true}}), C);
  if(D.comp&&D.comp.data&&D.comp.data.length)
    Plotly.newPlot('pt-comp', D.comp.data, merge(L,D.comp.layout), C);
  {clinical_js}
}})();
</script>
</body>
</html>"""


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

    # Generate dynamic insights last so they can read all of chart_data
    data["insights"] = generate_dynamic_insights(result, data)

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
        ("__INSIGHTS_JSON__",        json.dumps(chart_data.get("insights", {}), allow_nan=False)),
        ("__PLOTLY_SCRIPT__",        _PLOTLY_JS_INLINE),
    ]:
        html = html.replace(placeholder, value)
    return html
