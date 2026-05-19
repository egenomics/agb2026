"""charts/renderer.py — HTML report rendering (template filling and per-patient HTML).

This module is concerned only with presentation: it receives pre-computed chart
payloads and injects them into templates.  All analytical logic lives upstream
in orchestrator.py and insights.py.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .clinical import build_clinical_slope
from .config import BASE_CONFIG, BASE_LAYOUT, THEME
from .individual import build_faceted_composition, build_patient_radar_profiles

_PLOTLY_JS_PATH = Path(__file__).parent / "plotly.min.js"

# Loaded on first render call, not at import time.  This avoids downloading
# files or printing warnings when the module is merely imported (e.g. in tests
# that only exercise parsers) and prevents unexpected I/O on HPC nodes.
_PLOTLY_JS_INLINE: str | None = None
_TEMPLATE: str | None = None
_PATIENT_TEMPLATE: str | None = None


def _load_plotly_js() -> str:
    """Load Plotly.js inline for offline HPC use. Fail fast if the bundle is missing."""
    if not _PLOTLY_JS_PATH.exists():
        raise FileNotFoundError(
            f"Bundled Plotly.js not found: {_PLOTLY_JS_PATH}\n"
            "Commit plotly.min.js to git for offline HPC use:\n"
            "  curl -fsSL https://cdn.plot.ly/plotly-2.35.2.min.js \\\n"
            "       -o modules/groupD/microsee_report/report_generator/charts/plotly.min.js\n"
            "  git add modules/groupD/microsee_report/report_generator/charts/plotly.min.js"
        )
    return f"<script>{_PLOTLY_JS_PATH.read_text(encoding='utf-8')}</script>"


def _get_plotly_js() -> str:
    global _PLOTLY_JS_INLINE
    cached = _PLOTLY_JS_INLINE
    if cached is None:
        cached = _load_plotly_js()
        _PLOTLY_JS_INLINE = cached
    return cached


def _get_template() -> str:
    global _TEMPLATE
    cached = _TEMPLATE
    if cached is None:
        cached = (Path(__file__).parent / "template.html").read_text(encoding="utf-8")
        _TEMPLATE = cached
    return cached


def _get_patient_template() -> str:
    global _PATIENT_TEMPLATE
    cached = _PATIENT_TEMPLATE
    if cached is None:
        cached = (
            Path(__file__).parent / "patient_template.html"
        ).read_text(encoding="utf-8")
        _PATIENT_TEMPLATE = cached
    return cached


# ── Per-patient stability helper ──────────────────────────────────────────────

def _patient_bc(rows: list[dict[str, Any]], taxa: list[str], patient: str) -> float:
    """Bray-Curtis dissimilarity between T0 and T84 for one patient."""
    from .preprocessing import get_patient_timepoints
    r0, r84 = get_patient_timepoints(rows, patient)
    if r0 is None or r84 is None or not taxa:
        return 0.0
    ab0  = np.array([float(r0.get(t) or 0) for t in taxa])
    ab84 = np.array([float(r84.get(t) or 0) for t in taxa])
    denom = float(np.sum(ab0 + ab84))
    return round(float(np.sum(np.abs(ab0 - ab84))) / denom, 3) if denom > 0 else 0.0


# ── Per-patient personalised insight block ────────────────────────────────────

def _build_patient_insights_html(
    patient_id: str,
    p_rows: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    taxa: list[str],
    bc: float,
    pat_profile: dict[str, Any],
    has_clinical: bool,
) -> str:
    group   = str(p_rows[0].get("base_group") or p_rows[0].get("group") or "")
    t0_row  = next((r for r in p_rows if (r.get("time") or 0) == 0), None)
    t84_row = next((r for r in p_rows if (r.get("time") or 0) > 0), None)
    items: list[str] = []

    # 1 — Shannon diversity change
    if t0_row and t84_row:
        sh0  = float(t0_row.get("shannon") or 0)
        sh84 = float(t84_row.get("shannon") or 0)
        if sh0 > 0:
            delta_sh = sh84 - sh0
            pct_sh   = abs(delta_sh / sh0) * 100
            arrow    = "↑" if delta_sh > 0 else "↓"
            color    = "#D97A3A" if delta_sh > 0 else "#4A7ED4"
            verb     = "increased" if delta_sh > 0 else "decreased"
            items.append(
                f"<li><strong>Shannon diversity</strong> {verb} from "
                f"<strong>{sh0:.2f}</strong> → "
                f'<strong style="color:{color}">{sh84:.2f}</strong> '
                f'<span style="color:{color}">({arrow}{pct_sh:.0f}%)</span>.</li>'
            )

    # 2 — Largest taxon shifts (top 2)
    table = pat_profile.get("table", [])
    ranked = sorted(
        [r for r in table if r.get("delta") is not None],
        key=lambda r: abs(float(r["delta"])),
        reverse=True,
    )
    for entry in ranked[:2]:
        d   = float(entry["delta"])
        if abs(d) < 0.5:
            break
        fam   = entry["family"]
        dstr  = f"+{d:.1f}%" if d > 0 else f"{d:.1f}%"
        color = "#D97A3A" if d > 0 else "#4A7ED4"
        verb  = "increased" if d > 0 else "decreased"
        items.append(
            f"<li><strong>{fam}</strong> "
            f'<span style="color:{color}">{verb} by {dstr}</span> between visits.</li>'
        )

    # 3 — Stability rank within the same group
    group_patients = sorted(set(
        r["patient"] for r in rows
        if str(r.get("base_group") or r.get("group") or "") == group
    ))
    if len(group_patients) > 1:
        bc_scores = sorted(_patient_bc(rows, taxa, p) for p in group_patients)
        n_more_stable = sum(1 for v in bc_scores if v < bc)
        n_total       = len(bc_scores)
        group_mean    = round(float(np.mean(bc_scores)), 3)
        comparison    = (
            "more stable than" if bc < group_mean else
            "less stable than" if bc > group_mean else
            "as stable as"
        )
        items.append(
            f"<li>Microbiome stability score <strong>{bc:.3f}</strong> — "
            f"{comparison} the group average ({group_mean:.3f}). "
            f"More stable than {n_more_stable} of {n_total} patients in the {group} group.</li>"
        )

    # 4 — Clinical outcomes
    if has_clinical and t0_row and t84_row:
        mwt0  = float(t0_row.get("sixmwt") or 0)
        mwt84 = float(t84_row.get("sixmwt") or 0)
        if mwt0 > 0 and mwt84 > 0:
            dmwt  = mwt84 - mwt0
            color = "#D97A3A" if dmwt > 0 else "#4A7ED4"
            arrow = "↑" if dmwt > 0 else "↓"
            note  = " — physical function improved." if dmwt > 0 else " — physical function declined."
            items.append(
                f"<li><strong>6-Minute Walk Test</strong>: {mwt0:.0f} m → "
                f'<strong style="color:{color}">{mwt84:.0f} m</strong> '
                f'<span style="color:{color}">({arrow}{abs(dmwt):.0f} m)</span>{note}</li>'
            )
        il0  = float(t0_row.get("il18") or 0)
        il84 = float(t84_row.get("il18") or 0)
        if il0 > 0 and il84 > 0:
            dil   = il84 - il0
            color = "#4A7ED4" if dil < 0 else "#D97A3A"
            arrow = "↓" if dil < 0 else "↑"
            note  = " — inflammation marker reduced." if dil < 0 else " — inflammation marker increased."
            items.append(
                f"<li><strong>IL-18</strong>: {il0:.1f} → "
                f'<strong style="color:{color}">{il84:.1f} pg/mL</strong> '
                f'<span style="color:{color}">({arrow}{abs(dil):.1f})</span>{note}</li>'
            )

    if not items:
        return ""

    return (
        '<div class="insight">'
        f"<strong>Personalised findings — {patient_id} ({group}):</strong>"
        '<ul style="margin-top:10px;padding-left:18px;line-height:2.0">'
        + "".join(items)
        + "</ul></div>"
    )


# ── Per-patient report ────────────────────────────────────────────────────────

def render_patient_html(
    patient_id: str,
    result: Any,
    *,
    radar_profiles: dict[str, Any] | None = None,
) -> str:
    """Generate a self-contained per-patient HTML report.

    Contains: stability score, radar T0 vs T84 vs group mean, composition bars,
    clinical outcomes (if available), and a plain-language interpretation.

    Pass ``radar_profiles`` from :func:`build_patient_radar_profiles` when rendering
    many patients in one process (avoids recomputing all profiles per patient).
    """
    rows   = [r.model_dump() for r in result.rows]
    taxa   = result.taxa
    p_rows = [r for r in rows if r["patient"] == patient_id]
    if not p_rows:
        raise ValueError(f"Patient {patient_id!r} not found in result")

    bc    = _patient_bc(rows, taxa, patient_id)
    group: str = str(p_rows[0].get("base_group") or p_rows[0].get("group") or "")
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

    if radar_profiles is None:
        radar_profiles = build_patient_radar_profiles(rows, taxa)
    pat_profile = radar_profiles["profiles"].get(patient_id, {})
    patient_chart_data: dict[str, Any] = {
        "radar":       pat_profile.get("traces", []),
        "radar_table": pat_profile.get("table", []),
        "comp":        build_faceted_composition(p_rows, taxa),
        "bc":          bc,
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

    insights_html = _build_patient_insights_html(
        patient_id, p_rows, rows, taxa, bc, pat_profile, result.has_clinical,
    )

    replacements: list[tuple[str, str]] = [
        ("__PATIENT_ID__",        patient_id),
        ("__GROUP__",             group),
        ("__BC__",                f"{bc:.3f}"),
        ("__STAB_LABEL__",        stab_label),
        ("__BC_STABILITY_TEXT__", bc_stability_text),
        ("__CLINICAL_BLOCK__",    clinical_block),
        ("__CLINICAL_JS__",       clinical_js),
        ("__INSIGHTS_HTML__",     insights_html),
        ("__DATA_JSON__",         json.dumps(patient_chart_data, allow_nan=False)),
        ("__LAYOUT_JSON__",       json.dumps(BASE_LAYOUT, allow_nan=False)),
        ("__CONFIG_JSON__",       json.dumps(BASE_CONFIG,  allow_nan=False)),
        ("__FONT__",              THEME["font"]),
        ("__BG__",                THEME["bg"]),
        ("__TEXT__",              THEME["text"]),
        ("__PLOTLY_SCRIPT__",     _get_plotly_js()),
    ]
    html = _get_patient_template()
    for placeholder, value in replacements:
        html = html.replace(placeholder, value)
    return html


# ── Cohort report ─────────────────────────────────────────────────────────────

def render_html(chart_data: dict[str, Any]) -> str:
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

    html = _get_template()
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
        ("__PLOTLY_SCRIPT__",        _get_plotly_js()),
    ]:
        html = html.replace(placeholder, value)
    return html
