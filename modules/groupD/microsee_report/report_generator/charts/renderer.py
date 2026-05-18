"""charts/renderer.py — HTML report rendering (template filling and per-patient HTML).

This module is concerned only with presentation: it receives pre-computed chart
payloads and injects them into templates.  All analytical logic lives upstream
in orchestrator.py and insights.py.
"""
from __future__ import annotations
import json
import numpy as np
from pathlib import Path
from typing import Any

# Load Plotly.js from the bundled file so the report works offline on HPC nodes.
_PLOTLY_JS_PATH = Path(__file__).parent / "plotly.min.js"
if _PLOTLY_JS_PATH.exists():
    _PLOTLY_JS_INLINE = f"<script>{_PLOTLY_JS_PATH.read_text(encoding='utf-8')}</script>"
else:
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

_TEMPLATE         = (Path(__file__).parent / "template.html").read_text(encoding="utf-8")
_PATIENT_TEMPLATE = (Path(__file__).parent / "patient_template.html").read_text(encoding="utf-8")

from .config import BASE_LAYOUT, BASE_CONFIG, THEME
from .individual import build_patient_radar, build_faceted_composition
from .clinical   import build_clinical_slope
from .distances  import bray_curtis_matrix


# ── Per-patient stability helper ──────────────────────────────────────────────

def _patient_bc(rows: list[dict], taxa: list[str], patient: str) -> float:
    """Bray-Curtis dissimilarity between T0 and T84 for one patient."""
    from .preprocessing import get_patient_timepoints
    r0, r84 = get_patient_timepoints(rows, patient)
    if r0 is None or r84 is None or not taxa:
        return 0.0
    ab0  = np.array([float(r0.get(t) or 0) for t in taxa])
    ab84 = np.array([float(r84.get(t) or 0) for t in taxa])
    denom = float(np.sum(ab0 + ab84))
    return round(float(np.sum(np.abs(ab0 - ab84))) / denom, 3) if denom > 0 else 0.0


# ── Per-patient report ────────────────────────────────────────────────────────

def render_patient_html(patient_id: str, result: Any) -> str:
    """Generate a self-contained per-patient HTML report.

    Contains: stability score, radar T0 vs T84 vs group mean, composition bars,
    clinical outcomes (if available), and a plain-language interpretation.
    """
    rows   = [r.model_dump() for r in result.rows]
    taxa   = result.taxa
    p_rows = [r for r in rows if r["patient"] == patient_id]
    if not p_rows:
        raise ValueError(f"Patient {patient_id!r} not found in result")

    bc    = _patient_bc(rows, taxa, patient_id)
    group = p_rows[0].get("base_group", p_rows[0]["group"])
    stab_label = (
        "very stable"         if bc < 0.10 else
        "stable"              if bc < 0.20 else
        "moderately variable" if bc < 0.35 else
        "highly variable"
    )
    bc_stability_text = (
        "A low score indicates consistent microbial community structure."
        if bc < 0.2 else
        "This degree of change is common in dietary intervention studies."
    )

    patient_chart_data: dict[str, Any] = {
        "radar": build_patient_radar(rows, taxa),
        "comp":  build_faceted_composition(p_rows, taxa),
        "bc":    bc,
    }
    if result.has_clinical:
        patient_chart_data["sixmwt"] = build_clinical_slope(
            [r for r in rows if float(r.get("sixmwt") or 0) > 0], "sixmwt"
        )
        patient_chart_data["il18"] = build_clinical_slope(
            [r for r in rows if float(r.get("il18") or 0) > 0], "il18"
        )

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

    insight_text = (
        f"Patient {patient_id} ({group}) had a microbiome stability score of {bc:.3f}. "
        + ("Changes were minimal and within normal day-to-day variation." if bc < 0.2
           else "The largest compositional shift was detected, warranting closer monitoring.")
        + (" Clinical markers were available for this patient and are shown above."
           if result.has_clinical else "")
    )

    html = _PATIENT_TEMPLATE
    for placeholder, value in [
        ("__PATIENT_ID__",        patient_id),
        ("__GROUP__",             group),
        ("__BC__",                f"{bc:.3f}"),
        ("__STAB_LABEL__",        stab_label),
        ("__BC_STABILITY_TEXT__", bc_stability_text),
        ("__CLINICAL_BLOCK__",    clinical_block),
        ("__CLINICAL_JS__",       clinical_js),
        ("__INSIGHT_TEXT__",      insight_text),
        ("__DATA_JSON__",         json.dumps(patient_chart_data, allow_nan=False)),
        ("__LAYOUT_JSON__",       json.dumps(BASE_LAYOUT, allow_nan=False)),
        ("__CONFIG_JSON__",       json.dumps(BASE_CONFIG,  allow_nan=False)),
        ("__FONT__",              THEME["font"]),
        ("__BG__",                THEME["bg"]),
        ("__TEXT__",              THEME["text"]),
        ("__PLOTLY_SCRIPT__",     _PLOTLY_JS_INLINE),
    ]:
        html = html.replace(placeholder, value)
    return html


# ── Cohort report ─────────────────────────────────────────────────────────────

def render_html(chart_data: dict) -> str:
    """Render the main cohort HTML report by filling placeholders in template.html."""
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
