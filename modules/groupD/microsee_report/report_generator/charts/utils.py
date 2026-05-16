"""charts/utils.py — shared colour helpers and timepoint utilities."""

from __future__ import annotations
from .config import GROUP_COLORS, FALLBACK_COLOR, TAXA_COLORS, _EXTRA_TAXA_PALETTE

_taxon_color_cache: dict[str, str] = {}
_extra_index = 0


def group_color(group: str, all_groups: list[str]) -> str:
    try:
        return GROUP_COLORS[sorted(all_groups).index(group) % len(GROUP_COLORS)]
    except ValueError:
        return FALLBACK_COLOR


def taxon_color(taxon: str) -> str:
    global _extra_index
    if taxon in TAXA_COLORS:
        return TAXA_COLORS[taxon]
    if taxon not in _taxon_color_cache:
        _taxon_color_cache[taxon] = _EXTRA_TAXA_PALETTE[_extra_index % len(_EXTRA_TAXA_PALETTE)]
        _extra_index += 1
    return _taxon_color_cache[taxon]


def hex_rgba(hex_color: str, alpha: float) -> str:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def _base_group_color(base_group: str, base_groups: list[str]) -> str:
    """Color for a base group (EAA, Whey) using the same GROUP_COLORS palette."""
    try:
        return GROUP_COLORS[sorted(base_groups).index(base_group) % len(GROUP_COLORS)]
    except ValueError:
        return FALLBACK_COLOR


def _sorted_timepoints(rows: list[dict]) -> list[str]:
    """Return unique timepoints sorted by their numeric time value."""
    tp_time: dict[str, int] = {}
    for r in rows:
        tp = r.get("timepoint")
        if tp:
            tp_time[tp] = int(r.get("time") or 0)
    return sorted(tp_time.keys(), key=lambda tp: tp_time[tp])
