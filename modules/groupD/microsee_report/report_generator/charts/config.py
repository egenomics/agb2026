"""charts/config.py — colour constants and Plotly layout/config defaults."""

from __future__ import annotations

from typing import Any

TAXA_COLORS: dict[str, str] = {
    "Bacteroidaceae": "#4A7ED4",
    "Lachnospiraceae": "#D84E6A",
    "Ruminococcaceae": "#2FA896",
    "Prevotellaceae": "#D97A3A",
    "Rikenellaceae": "#9058C4",
    "Enterobacteriaceae": "#C4960A",
    "Oscillospiraceae": "#8CAD70",
    "Tannerellaceae": "#E87AA0",
    "Akkermansiaceae": "#5B8C5A",
}

# Extended palette for taxa families not listed above — cycles instead of all-grey.
EXTRA_TAXA_PALETTE = [
    "#6CB4E4",
    "#F4A261",
    "#E76F51",
    "#264653",
    "#2A9D8F",
    "#E9C46A",
    "#A8DADC",
    "#457B9D",
    "#1D3557",
    "#F1FAEE",
    "#B5838D",
    "#6D6875",
    "#EDC4B3",
    "#E8A598",
    "#CFBAF0",
    "#A3C4F3",
    "#90DBF4",
    "#8EECF5",
    "#98F5E1",
    "#B9FBC0",
]

GROUP_COLORS = ["#D97A3A", "#C4960A", "#4A7ED4", "#89BAF0", "#2FA896", "#9058C4"]
FALLBACK_COLOR = "#aaaaaa"

THEME = {
    "bg": "#FEF3EC",
    "paper": "#FFFFFF",
    "grid": "rgba(196,160,140,0.12)",
    "text": "#6B3A2A",
    "text2": "#8B5860",
    "font": "Nunito, system-ui, sans-serif",
}

BASE_LAYOUT: dict[str, Any] = {
    "autosize": True,
    "paper_bgcolor": THEME["paper"],
    "plot_bgcolor": THEME["bg"],
    "font": {"family": THEME["font"], "color": THEME["text"], "size": 11},
    "margin": {"l": 52, "r": 16, "t": 32, "b": 40},
    "hoverlabel": {"bgcolor": THEME["paper"], "font": {"size": 11, "family": THEME["font"]}},
    "xaxis": {
        "gridcolor": THEME["grid"],
        "linecolor": THEME["grid"],
        "zerolinecolor": THEME["grid"],
        "tickfont": {"size": 10, "color": THEME["text2"]},
    },
    "yaxis": {
        "gridcolor": THEME["grid"],
        "linecolor": THEME["grid"],
        "zerolinecolor": THEME["grid"],
        "tickfont": {"size": 10, "color": THEME["text2"]},
    },
    "legend": {
        "font": {"size": 10},
        "bgcolor": "rgba(0,0,0,0)",
        "orientation": "h",
        "yanchor": "bottom",
        "y": -0.28,
        "xanchor": "left",
        "x": 0,
    },
}

BASE_CONFIG: dict[str, Any] = {
    "displaylogo": False,
    "responsive": True,
    "modeBarButtonsToRemove": ["lasso2d", "select2d", "autoScale2d"],
    "toImageButtonOptions": {"format": "png", "scale": 2},
}
