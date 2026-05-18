"""charts — public API for the MicroSee report generator chart subpackage.

compute_chart_data  — orchestrates all chart-data computation (orchestrator.py)
render_html         — fills the HTML template with chart payloads (renderer.py)

Both names are re-exported here for backward compatibility: callers that do
    from report_generator.charts import compute_chart_data, render_html
continue to work without any changes.
"""
from .orchestrator import compute_chart_data, ReportConfig  # noqa: F401
from .renderer     import render_html                        # noqa: F401

__all__ = ["compute_chart_data", "render_html", "ReportConfig"]
