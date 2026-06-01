#!/usr/bin/env python3
"""
MicroSee Album Report Generator
================================
Folder-driven report builder for the Group D microbiome visualisation module.

Scans an input directory of PNGs (staged by Nextflow channels) and assembles
them into self-contained HTML reports matching the MicroSee visual identity.

Output:
  - index.html              main report: ice-cream splash, exploratory cohort
                            plots, and a grid of patient cards linking out to
                            individual patient pages.
  - patients/<ID>.html      one standalone, self-contained, separately
                            downloadable report per patient.

Design philosophy (inherited from MicroSee):
  - GRACEFUL SKIPPING: only renders PNGs that exist. Missing plots do not
    appear. No crash, no guessing filenames.
  - SELF-CONTAINED: every image is base64-embedded. Each HTML opens anywhere
    offline with no sibling files (cluster / HPC safe). Each patient page is
    independently downloadable.
  - FOLDER-DRIVEN: no hardcoded filenames; whatever PNGs land in the folders
    get displayed.

Expected input layout (any subset may be present):
    <input>/exploratory/*.png
    <input>/explanatory/alpha_patient/*.png
    <input>/explanatory/PCoA_patient/*.png
    <input>/explanatory/<anything>/*.png

Usage:
    python build_album.py --input build_input --outdir report
    # then open report/index.html
"""

from __future__ import annotations

import argparse
import base64
import csv
import html
import re
import sys
from datetime import datetime
from pathlib import Path

BG = "#FEF3EC"
PAPER = "#FFFFFF"
TEXT = "#6B3A2A"
TEXT2 = "#8B5860"
ACCENT = "#D97A3A"
FONT = "Nunito, system-ui, sans-serif"

PATIENT_RE = re.compile(r"(patient\d+|ERR\d{6,})", re.IGNORECASE)

# Friendly labels for known plot filenames (matched on the file stem).
# Anything not listed falls back to a prettified version of the filename,
# so unknown/new plots still render with a readable title.
PRETTY = {
    # exploratory / cohort
    "Table1_Demographics_Plot": "Cohort Demographics",
    "demographic_table": "Cohort Demographics",
    "pca_samples_bacteria": "PCA — Samples vs Bacteria",
    "pca_biplot_healthy_vs_disease": "PCA Biplot — Healthy vs Disease",
    "pca_individuals": "PCA — Individuals",
    "pca_variables": "PCA — Variables",
    "pca_scree_plot": "PCA — Scree Plot",
    "pca_pc1_boxplot": "PC1 Boxplot",
    "pca_pc1_density": "PC1 Density",
    "parallel_plot_relative_abundance": "Relative Abundance (Parallel)",
    "parallel_plot_relative_abundance_healthygrouped": "Relative Abundance (Parallel)",
    "heatmap_samples_bacteria": "Heatmap — Samples vs Bacteria",
    "volcano_plot_ancombc2": "Volcano Plot (ANCOM-BC2)",
    "volcano_plot": "Volcano Plot",
    "beta_diversity_tree": "Beta Diversity Dendrogram",
    "clinical_cooccurrence_heatmap": "Clinical Co-occurrence Heatmap",
    "clinical_association_map": "Clinical Association Map",
}

# Preferred display order for cohort plots (unknowns appended after, A–Z).
EXPLORATORY_ORDER = [
    "demographic_table", "Table1_Demographics_Plot",
    "pca_biplot_healthy_vs_disease", "pca_samples_bacteria",
    "pca_individuals", "pca_variables", "pca_scree_plot",
    "pca_pc1_boxplot", "pca_pc1_density", "PC1_bacterias", "boxplot_PC1",
    "parallel_plot_relative_abundance_healthygrouped",
    "parallel_plot_relative_abundance",
    "heatmap_samples_bacteria",
    "volcano_plot_ancombc2", "volcano_plot",
    "beta_diversity_tree", "clinical_cooccurrence_heatmap",
    "clinical_association_map",
]

# Per-patient plot-type label keyed by the directory the PNG was found under.
PER_PATIENT_GROUP_LABEL = {
    "alpha_patient": "Alpha Diversity Distribution",
    "PCoA_patient": "PCoA — Patient Highlighted",
    "overview_sample": "Patient Overview",
    "Virulence_analysis": "Virulence Profile",
}


def b64_img(path: Path) -> str:
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{data}"


def pretty_name(stem: str) -> str:
    if stem in PRETTY:
        return PRETTY[stem]
    cleaned = re.sub(r"[_\-]+", " ", stem).strip()
    return cleaned[:1].upper() + cleaned[1:] if cleaned else stem


def _per_patient_label(png: Path) -> str:
    """Choose a label for a per-patient plot based on the directory it sits in."""
    for part in reversed(png.parts):
        if part in PER_PATIENT_GROUP_LABEL:
            return PER_PATIENT_GROUP_LABEL[part]
    # fall back to parent folder name, prettified
    parent = png.parent.name
    if parent in PER_PATIENT_GROUP_LABEL:
        return PER_PATIENT_GROUP_LABEL[parent]
    return pretty_name(parent if parent else png.stem)


def _iter_pngs(root: Path):
    """Yield every .png under root, recursively. Non-PNGs (.rds, .html, .qzv) ignored."""
    if not root.is_dir():
        return
    yield from sorted(root.rglob("*.png"))


def collect_patient_tsvs(root: Path) -> dict:
    """
    Per-patient TSV tables, grouped by patient id. Scans recursively for .tsv
    files whose filename contains a patient id, parses them with csv, and groups
    by pid. Returns {pid: {label: [[header_row], [data_row], ...]}}.
    Non-TSV files and files without a patient id in the name are ignored.
    Malformed or empty files are skipped without crashing.
    """
    patients: dict = {}
    seen: set = set()
    for tsv in sorted(root.rglob("*.tsv")):
        rp = tsv.resolve()
        if rp in seen:
            continue
        m = PATIENT_RE.search(tsv.name)
        if not m:
            continue
        seen.add(rp)
        pid = m.group(1)
        stem_clean = re.sub(PATIENT_RE, "", tsv.stem).strip("_-")
        label = pretty_name(stem_clean) if stem_clean else pretty_name(tsv.stem)
        try:
            with tsv.open(newline="", encoding="utf-8") as fh:
                rows = list(csv.reader(fh, delimiter="\t"))
            if not rows:
                continue
        except Exception:
            continue
        bucket = patients.setdefault(pid, {})
        key = label
        n = 2
        while key in bucket:
            key = f"{label} ({n})"
            n += 1
        bucket[key] = rows
    return patients


def tsv_card(title: str, rows: list) -> str:
    """Render a parsed TSV (list of row-lists) as a styled HTML table card."""
    if not rows:
        return ""
    header, *body = rows
    th_cells = "".join(
        f'<th style="text-align:left;padding:6px 8px;'
        f'border-bottom:1px solid rgba(196,160,140,.2);'
        f'font-size:11px;color:{TEXT2};font-weight:600;white-space:nowrap">{html.escape(h)}</th>'
        for h in header
    )
    td_rows = "".join(
        "<tr>" + "".join(
            f'<td style="padding:5px 8px;border-bottom:1px solid rgba(196,160,140,.1);'
            f'font-size:11px;color:{TEXT};vertical-align:top">{html.escape(cell)}</td>'
            for cell in row
        ) + "</tr>"
        for row in body if row
    )
    return (
        f'<div class="chart-card">'
        f'<div class="chart-title">{title}</div>'
        f'<div style="overflow-x:auto">'
        f'<table style="width:100%;border-collapse:collapse">'
        f'<thead><tr style="background:{BG}">{th_cells}</tr></thead>'
        f'<tbody>{td_rows}</tbody>'
        f'</table></div></div>'
    )


def collect_cohort_tsvs(root: Path) -> list:
    """
    Cohort (non per-patient) TSV tables. Scans recursively for .tsv files whose
    filename does NOT contain a patient id and that are not under a per-patient
    subtree. Mirrors the exclusion logic of collect_exploratory.
    Returns a list of (label, rows) pairs sorted by filename.
    Malformed or empty files are skipped silently.
    """
    explanatory_markers = {"explanatory", "explanatory_report",
                           "PCoA_patient", "alpha_patient", "overview_sample",
                           "Virulence_analysis"}

    def under_explanatory(p: Path) -> bool:
        return any(part in explanatory_markers for part in p.parts)

    results = []
    seen: set = set()
    for tsv in sorted(root.rglob("*.tsv")):
        rp = tsv.resolve()
        if rp in seen:
            continue
        if PATIENT_RE.search(tsv.name):
            continue
        if under_explanatory(tsv):
            continue
        seen.add(rp)
        try:
            with tsv.open(newline="", encoding="utf-8") as fh:
                rows = list(csv.reader(fh, delimiter="\t"))
            if not rows:
                continue
        except Exception:
            continue
        results.append((pretty_name(tsv.stem), rows))
    return results


def collect_exploratory(root: Path) -> list[Path]:
    """
    Cohort (non per-patient) PNGs. Searches the whole input tree recursively
    and keeps any PNG whose filename does NOT contain a patient (ERR) id, and
    that does not live under an explanatory / per-patient subtree.
    """
    pngs: list[Path] = []
    seen = set()

    # INSULATION: Define exclusions using lowercase file stems to match sorting mechanics
    excluded_stems = {
        "pca_individuals", 
        "pca_scree_plot"
    }

    # Directories whose contents are per-patient and must be excluded here.
    explanatory_markers = {"explanatory", "explanatory_report",
                           "PCoA_patient", "alpha_patient", "overview_sample",
                           "Virulence_analysis"}

    def under_explanatory(p: Path) -> bool:
        return any(part in explanatory_markers for part in p.parts)

    for p in _iter_pngs(root):
        rp = p.resolve()
        if rp in seen:
            continue
        if PATIENT_RE.search(p.name):      # per-patient → not cohort
            continue
        if under_explanatory(p):           # inside a per-patient subtree → skip
            continue
        
        # FIX: Check against a case-insensitive conversion of the file stem
        if p.stem.lower() in excluded_stems:
            continue
            
        seen.add(rp)
        pngs.append(p)

    rank = {name: i for i, name in enumerate(EXPLORATORY_ORDER)}
    return sorted(pngs, key=lambda p: (rank.get(p.stem, 999), p.stem))

def collect_explanatory(root: Path) -> dict:
    """
    Per-patient PNGs, grouped by patient id. Any PNG whose filename contains an
    ERR id is treated as per-patient, regardless of how deeply it is nested
    (PCoA_patient/, alpha_patient/, overview_sample/, explanatory_report/...).
    """
    patients: dict = {}
    seen = set()
    for png in _iter_pngs(root):
        if png.resolve() in seen:
            continue
        m = PATIENT_RE.search(png.name)
        if not m:
            continue
        seen.add(png.resolve())
        pid = m.group(1)
        label = _per_patient_label(png)
        # if two plots map to the same label for a patient, disambiguate
        bucket = patients.setdefault(pid, {})
        key = label
        n = 2
        while key in bucket:
            key = f"{label} ({n})"
            n += 1
        bucket[key] = png
    return dict(sorted(patients.items()))


def css() -> str:
    return f"""
*{{box-sizing:border-box;margin:0;padding:0}}
html{{scroll-behavior:smooth}}
body{{font-family:{FONT};background:{BG};color:{TEXT};display:flex;min-height:100vh}}
#sidebar{{width:210px;min-height:100vh;background:{PAPER};border-right:1px solid rgba(196,160,140,.2);
  padding:24px 0;position:sticky;top:0;height:100vh;overflow-y:auto;flex-shrink:0}}
#sidebar .logo{{padding:0 20px 20px;font-size:18px;font-weight:700;color:{TEXT};
  border-bottom:1px solid rgba(196,160,140,.15);margin-bottom:12px}}
#sidebar .logo span{{font-size:11px;font-weight:400;color:{TEXT2};display:block;margin-top:2px}}
#sidebar a{{display:block;padding:8px 20px;font-size:13px;color:{TEXT2};text-decoration:none;
  border-left:3px solid transparent;transition:all .15s}}
#sidebar a:hover,#sidebar a.active{{color:{TEXT};border-left-color:{ACCENT};background:rgba(217,122,58,.06)}}
.nav-group-label{{font-size:10px;color:{TEXT2};text-transform:uppercase;letter-spacing:.07em;
  font-weight:600;padding:12px 20px 4px;border-top:1px solid rgba(196,160,140,.15);margin-top:4px}}
#wrap{{flex:1;display:flex;flex-direction:column;overflow-x:hidden}}
#topbar{{background:{PAPER};border-bottom:1px solid rgba(196,160,140,.2);
  padding:12px 24px;display:flex;gap:32px;align-items:center;position:sticky;top:0;z-index:50}}
.stat{{display:flex;flex-direction:column}}
.stat-val{{font-size:20px;font-weight:700;color:{TEXT}}}
.stat-lbl{{font-size:10px;color:{TEXT2};text-transform:uppercase;letter-spacing:.05em}}
.topbtn{{font-size:12px;font-weight:700;padding:8px 16px;border:none;border-radius:20px;
  background:{ACCENT};color:#fff;cursor:pointer;font-family:{FONT};transition:all .12s;text-decoration:none;
  display:inline-flex;align-items:center}}
.topbtn:hover{{background:#c06a2e}}
.topbtn.ghost{{background:transparent;color:{TEXT2};border:1px solid rgba(196,160,140,.4)}}
.topbtn.ghost:hover{{background:rgba(217,122,58,.08);color:{TEXT}}}
.ml-auto{{margin-left:auto}}
.btn-row{{display:flex;gap:8px;align-items:center}}
#main{{flex:1;padding:24px;overflow-x:hidden}}
section{{margin-bottom:40px;scroll-margin-top:70px}}
.sec-header{{font-size:16px;font-weight:700;color:{TEXT};margin-bottom:10px;
  padding-bottom:8px;border-bottom:2px solid rgba(196,160,140,.2)}}
.sec-insight{{font-size:12px;color:{TEXT2};line-height:1.55;
  border-left:3px solid {ACCENT};padding:4px 0 4px 10px;margin-bottom:14px}}
.grid2{{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin-bottom:16px}}
@media(max-width:900px){{.grid2{{grid-template-columns:1fr}}}}
.chart-card{{background:{PAPER};border-radius:10px;padding:16px;box-shadow:0 1px 4px rgba(107,58,42,.08)}}
.chart-title{{font-size:13px;font-weight:700;color:{TEXT};margin-bottom:10px}}
.album-img{{width:100%;height:auto;border-radius:6px;cursor:zoom-in;display:block;background:#faf5f0}}
.pgrid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:14px}}
.pcard{{background:{PAPER};border-radius:10px;padding:18px 14px;text-decoration:none;
  box-shadow:0 1px 4px rgba(107,58,42,.08);transition:all .14s;display:flex;
  flex-direction:column;gap:8px;border:1px solid transparent}}
.pcard:hover{{transform:translateY(-2px);box-shadow:0 4px 14px rgba(107,58,42,.16);border-color:{ACCENT}}}
.pcard .pid{{font-size:14px;font-weight:700;color:{TEXT}}}
.pcard .pmeta{{font-size:11px;color:{TEXT2}}}
.pcard .open{{font-size:11px;font-weight:700;color:{ACCENT};margin-top:auto}}
#lb{{position:fixed;inset:0;background:rgba(62,26,14,.85);z-index:9999;display:none;
  align-items:center;justify-content:center;cursor:zoom-out;padding:30px}}
#lb.open{{display:flex}}
#lb-card{{background:{PAPER};border-radius:14px;padding:16px;cursor:default;
  max-width:92vw;max-height:92vh;display:flex;flex-direction:column;gap:12px;
  box-shadow:0 12px 50px rgba(0,0,0,.45);animation:cardPop .22s cubic-bezier(.34,1.56,.64,1)}}
@keyframes cardPop{{0%{{opacity:0;transform:scale(.9)}}100%{{opacity:1;transform:scale(1)}}}}
#lb-card img{{max-width:88vw;max-height:74vh;border-radius:8px;display:block;background:#faf5f0}}
.lb-head{{display:flex;align-items:center;gap:16px}}
#lb .cap{{color:{TEXT};font-size:14px;font-weight:700;flex:1}}
.lb-actions{{display:flex;gap:8px;align-items:center}}
.lb-btn{{font-size:12px;font-weight:700;padding:7px 14px;border:none;border-radius:18px;
  background:{ACCENT};color:#fff;cursor:pointer;font-family:{FONT};text-decoration:none;
  display:inline-flex;align-items:center;transition:all .12s}}
.lb-btn:hover{{background:#c06a2e}}
.lb-btn.ghost{{background:transparent;color:{TEXT2};border:1px solid rgba(196,160,140,.4)}}
.lb-btn.ghost:hover{{background:rgba(217,122,58,.08);color:{TEXT}}}
#lb .x{{position:absolute;top:18px;right:24px;color:#fff;font-size:30px;cursor:pointer;line-height:1;font-weight:300}}
@keyframes scoopPop{{0%{{opacity:0;transform:translateY(14px) scale(.75)}}
  65%{{transform:translateY(-2px) scale(1.05)}}100%{{opacity:1;transform:translateY(0) scale(1)}}}}
@keyframes fadeInUp{{from{{opacity:0;transform:translateY(8px)}}to{{opacity:1;transform:translateY(0)}}}}
@keyframes splashExit{{0%{{opacity:1;visibility:visible}}99%{{opacity:0;visibility:visible}}
  100%{{opacity:0;visibility:hidden}}}}
#splash{{position:fixed;top:0;left:0;width:100%;height:100%;background:#FEF3EC;z-index:10000;
  display:flex;flex-direction:column;align-items:center;justify-content:center;gap:6px;
  animation:splashExit .65s ease-in 3.2s forwards}}
#splash .sc{{opacity:0;transform-box:fill-box;transform-origin:center bottom}}
#splash .sc1{{animation:scoopPop .4s cubic-bezier(.34,1.56,.64,1) .35s forwards}}
#splash .sc2{{animation:scoopPop .4s cubic-bezier(.34,1.56,.64,1) .68s forwards}}
#splash .sc3{{animation:scoopPop .4s cubic-bezier(.34,1.56,.64,1) 1.01s forwards}}
#splash .sc4{{animation:scoopPop .4s cubic-bezier(.34,1.56,.64,1) 1.34s forwards}}
.splash-brand{{opacity:0;font-size:30px;font-weight:800;color:#3E1A0E;letter-spacing:-.3px;
  font-family:{FONT};animation:fadeInUp .55s ease-out 1.9s forwards}}
.splash-sub{{opacity:0;font-size:12px;color:#8B5860;letter-spacing:.04em;
  font-family:{FONT};animation:fadeInUp .55s ease-out 2.35s forwards}}
@media print{{
  #sidebar,#topbar,.topbtn,#splash{{display:none!important}}
  body{{display:block}} #main{{padding:0}}
  .grid2{{grid-template-columns:1fr 1fr}}
  .chart-card{{box-shadow:none;border:1px solid #eee;break-inside:avoid;page-break-inside:avoid}}
  .album-img{{cursor:default}}
}}
"""


def splash() -> str:
    return """
<div id="splash" role="presentation" aria-hidden="true">
  <svg width="110" height="200" viewBox="0 0 110 200" xmlns="http://www.w3.org/2000/svg">
    <polygon points="55,200 8,118 102,118" fill="#D97A3A"/>
    <line x1="55" y1="200" x2="31" y2="135" stroke="rgba(0,0,0,0.13)" stroke-width="1.2"/>
    <line x1="55" y1="200" x2="79" y2="135" stroke="rgba(0,0,0,0.13)" stroke-width="1.2"/>
    <line x1="18" y1="155" x2="92" y2="155" stroke="rgba(0,0,0,0.10)" stroke-width="1"/>
    <ellipse cx="55" cy="118" rx="47" ry="8" fill="#B85A18"/>
    <circle class="sc sc1" cx="55" cy="103" r="32" fill="#C4960A"/>
    <circle class="sc sc2" cx="55" cy="73"  r="28" fill="#2FA896"/>
    <circle class="sc sc3" cx="55" cy="47"  r="24" fill="#D84E6A"/>
    <circle class="sc sc4" cx="55" cy="25"  r="20" fill="#4A7ED4"/>
    <circle class="sc sc4" cx="62" cy="6" r="7" fill="#9058C4"/>
    <line x1="55" y1="6" x2="55" y2="14" stroke="#5A3010" stroke-width="1.5" stroke-linecap="round"/>
  </svg>
  <div class="splash-brand">MicroSee</div>
  <div class="splash-sub">Microbiome Visualisation Report</div>
</div>"""


def lightbox_and_js() -> str:
    return """
<div id="lb" onclick="closeLb(event)">
  <span class="x" onclick="closeLb(event)">&times;</span>
  <div id="lb-card" onclick="event.stopPropagation()">
    <div class="lb-head">
      <div class="cap" id="lb-cap"></div>
      <div class="lb-actions">
        <a id="lb-dl" class="lb-btn" download="">⬇ Download PNG</a>
        <button class="lb-btn ghost" onclick="closeLb(event)">Close</button>
      </div>
    </div>
    <img id="lb-img" src="" alt="">
  </div>
</div>
<script>
function zoom(src,cap,fname){
  document.getElementById('lb-img').src=src;
  document.getElementById('lb-cap').textContent=cap;
  var dl=document.getElementById('lb-dl');
  dl.href=src; dl.setAttribute('download',(fname||cap||'plot')+'.png');
  document.getElementById('lb').classList.add('open');
}
function closeLb(e){document.getElementById('lb').classList.remove('open');}
document.addEventListener('keydown',function(e){if(e.key==='Escape')closeLb(e);});
</script>"""


def img_card(title: str, src: str, fname: str = "") -> str:
    safe = title.replace("'", "")
    fsafe = (fname or title).replace("'", "")
    return f"""
      <div class="chart-card">
        <div class="chart-title">{title}</div>
        <img class="album-img" src="{src}" alt="{title}" loading="lazy"
             onclick="zoom(this.src,'{safe}','{fsafe}')">
      </div>"""


def build_patient_page(pid: str, groups: dict, tsv_groups: dict | None = None, standalone: bool = False) -> str:
    if tsv_groups is None:
        tsv_groups = {}
    png_cards = "\n".join(img_card(label, b64_img(path), path.stem) for label, path in groups.items())
    # Metadata table is large (key/value, many rows) → render full-width below the
    # grid. All other tables (e.g. Top genus) stay in the grid with the plots.
    grid_tables = {l: r for l, r in tsv_groups.items() if "metadata" not in l.lower()}
    wide_tables = {l: r for l, r in tsv_groups.items() if "metadata" in l.lower()}
    table_cards = "\n".join(tsv_card(label, rows) for label, rows in grid_tables.items())
    all_cards = png_cards + ("\n" + table_cards if table_cards else "")
    wide_cards = "\n".join(tsv_card(label, rows) for label, rows in wide_tables.items())
    n_items = len(groups) + len(tsv_groups)
    generated = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M")
    back_sidebar = "" if standalone else '<a href="../index.html">← Back to overview</a>'
    back_top = "" if standalone else '<a class="topbtn ghost" href="../index.html">← Overview</a>'
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>MicroSee — Patient {pid}</title>
<style>{css()}</style></head>
<body>
<nav id="sidebar">
  <div class="logo">MicroSee<span>Patient Report</span></div>
  {back_sidebar}
</nav>
<div id="wrap">
  <div id="topbar">
    <div class="stat"><div class="stat-val">{pid}</div><div class="stat-lbl">Patient</div></div>
    <div class="stat"><div class="stat-val">{n_items}</div><div class="stat-lbl">Items</div></div>
    <div class="btn-row ml-auto">
      {back_top}
      <button class="topbtn" onclick="window.print()">⬇ Download PDF</button>
    </div>
  </div>
  <div id="main">
    <section>
      <div class="sec-header">Patient {pid} — Individual Report</div>
      <div class="sec-insight">Per-patient diversity plots and tables. Click any plot to enlarge. Use Download PDF to save this patient's report on its own.</div>
      <div class="grid2">{all_cards}</div>
      {wide_cards}
    </section>
    <p style="font-size:11px;color:{TEXT2}">Generated {generated} · AGB 2026 - Group D</p>
  </div>
</div>
{lightbox_and_js()}
</body></html>"""


def build_index(exploratory: list, explanatory: dict, patient_tsvs: dict | None = None,
                cohort_tsvs: list | None = None) -> str:
    if patient_tsvs is None:
        patient_tsvs = {}
    if cohort_tsvs is None:
        cohort_tsvs = []
    n_explor = len(exploratory)
    n_patients = len(explanatory)
    n_cohort_tsvs = len(cohort_tsvs)
    n_total = n_explor + n_cohort_tsvs + sum(
        len(v) + len(patient_tsvs.get(pid, {})) for pid, v in explanatory.items()
    )

    nav_links = ['<a href="#exploratory" class="active">Exploratory</a>',
                 '<a href="#patients">Patient reports</a>']
    if explanatory:
        nav_links.append('<div class="nav-group-label">Patients</div>')
        for pid in explanatory:
            nav_links.append(f'<a href="patients/{pid}.html" target="_blank" rel="noopener">{pid}</a>')
    nav_html = "\n".join(nav_links)

    if exploratory or cohort_tsvs:
        png_cards = "\n".join(img_card(pretty_name(p.stem), b64_img(p), p.stem) for p in exploratory)
        tbl_cards = "\n".join(tsv_card(label, rows) for label, rows in cohort_tsvs)
        all_cohort = "\n".join(filter(None, [png_cards, tbl_cards]))
        explor_html = f"""
    <section id="exploratory">
      <div class="sec-header">Exploratory Report — Whole Cohort</div>
      <div class="sec-insight">Cohort-level overview across all samples. Click any plot to enlarge.</div>
      <div class="grid2">{all_cohort}</div>
    </section>"""
    else:
        explor_html = """
    <section id="exploratory">
      <div class="sec-header">Exploratory Report — Whole Cohort</div>
      <div class="sec-insight">No cohort-level plots found yet.</div>
    </section>"""

    pcards = []
    for pid, groups in explanatory.items():
        pcards.append(f"""
        <a class="pcard" href="patients/{pid}.html" target="_blank" rel="noopener">
          <div class="pid">{pid}</div>
          <div class="pmeta">{len(groups) + len(patient_tsvs.get(pid, {}))} item(s)</div>
          <div class="open">Open report ↗</div>
        </a>""")
    patients_html = f"""
    <section id="patients">
      <div class="sec-header">Individual Patient Reports</div>
      <div class="sec-insight">Each patient has a standalone report that opens in a new tab and can be downloaded separately as PDF.</div>
      <div class="pgrid">{''.join(pcards)}</div>
    </section>""" if explanatory else ""

    generated = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M")

    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>MicroSee Report</title>
<style>{css()}</style></head>
<body>
{splash()}
<nav id="sidebar">
  <div class="logo">MicroSee<span>Group D — Visualisation</span></div>
  {nav_html}
</nav>
<div id="wrap">
  <div id="topbar">
    <div class="stat"><div class="stat-val">{n_total}</div><div class="stat-lbl">Plots</div></div>
    <div class="stat"><div class="stat-val">{n_explor}</div><div class="stat-lbl">Exploratory</div></div>
    <div class="stat"><div class="stat-val">{n_patients}</div><div class="stat-lbl">Patients</div></div>
    <button class="topbtn ml-auto" onclick="window.print()">⬇ Download PDF</button>
  </div>
  <div id="main">
    {explor_html}
    {patients_html}
    <p style="font-size:11px;color:{TEXT2};margin-top:24px">Generated {generated} · AGB 2026 - Group D</p>
  </div>
</div>
{lightbox_and_js()}
<script>
const links=[...document.querySelectorAll('#sidebar a[href^="#"]')];
const secs=links.map(a=>document.querySelector(a.getAttribute('href'))).filter(Boolean);
window.addEventListener('scroll',()=>{{let cur=secs[0];
  for(const s of secs){{if(s.getBoundingClientRect().top<120)cur=s;}}
  links.forEach(a=>a.classList.toggle('active',a.getAttribute('href')==='#'+(cur&&cur.id)));}});
</script>
</body></html>"""


def main() -> int:
    ap = argparse.ArgumentParser(
        description="MicroSee album report generator",
        epilog="Modes: full report (default) | --patient ID (one patient only) "
               "| --explanatory-only (all patient pages, no cohort index)")
    ap.add_argument("--input", required=True,
                    help="Input dir with exploratory/ and/or explanatory/")
    ap.add_argument("--outdir", default="report", help="Output directory")
    ap.add_argument("--patient", default=None,
                    help="Generate ONLY this patient's standalone report (e.g. ERR1080200). "
                         "Output is a single self-contained HTML.")
    ap.add_argument("--explanatory-only", action="store_true",
                    help="Generate only the per-patient reports, skip the cohort index.")
    args = ap.parse_args()

    root = Path(args.input)
    if not root.is_dir():
        print(f"ERROR: input dir not found: {root}", file=sys.stderr)
        return 1

    explanatory = collect_explanatory(root)
    patient_tsvs = collect_patient_tsvs(root)
    cohort_tsvs = collect_cohort_tsvs(root)

    # ── Mode 1: single patient standalone ───────────────────────────────────
    if args.patient:
        groups = explanatory.get(args.patient)
        if not groups:
            print(f"ERROR: no plots found for patient {args.patient} under {root}/explanatory/",
                  file=sys.stderr)
            return 1
        out = Path(args.outdir)
        if out.suffix.lower() == ".html":
            out.parent.mkdir(parents=True, exist_ok=True)
            target = out
        else:
            out.mkdir(parents=True, exist_ok=True)
            target = out / f"{args.patient}.html"
        target.write_text(build_patient_page(args.patient, groups,
                                             tsv_groups=patient_tsvs.get(args.patient, {}),
                                             standalone=True),
                          encoding="utf-8")
        print(f"Wrote standalone patient report: {target}")
        return 0

    # ── Mode 2/3: patient pages (+ optional cohort index) ───────────────────
    outdir = Path(args.outdir)
    (outdir / "patients").mkdir(parents=True, exist_ok=True)

    for pid, groups in explanatory.items():
        (outdir / "patients" / f"{pid}.html").write_text(
            build_patient_page(pid, groups, tsv_groups=patient_tsvs.get(pid, {})),
            encoding="utf-8")
    print(f"Wrote {len(explanatory)} patient page(s) to {outdir}/patients/")

    if not args.explanatory_only:
        exploratory = collect_exploratory(root)
        (outdir / "index.html").write_text(
            build_index(exploratory, explanatory, patient_tsvs, cohort_tsvs), encoding="utf-8")
        print(f"Wrote {outdir}/index.html")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
