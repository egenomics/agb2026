"""charts — public API for the MicroSee report generator chart subpackage.

compute_chart_data   — orchestrates all chart-data computation (orchestrator.py)
ReportConfig         — section-selection config passed to compute_chart_data
render_html          — fills the cohort HTML template with chart payloads (renderer.py)
render_patient_html  — generates a self-contained per-patient HTML report (renderer.py)
"""

from .orchestrator import ReportConfig, compute_chart_data
from .renderer import render_html, render_patient_html

__all__ = ["ReportConfig", "compute_chart_data", "render_html", "render_patient_html"]
