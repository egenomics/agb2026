"""charts/utils.py — shared colour helpers and styling utilities."""
from __future__ import annotations
from .config import GROUP_COLORS, FALLBACK_COLOR, TAXA_COLORS, _EXTRA_TAXA_PALETTE


def group_color(group: str, all_groups: list[str]) -> str:
    try:
        return GROUP_COLORS[sorted(all_groups).index(group) % len(GROUP_COLORS)]
    except ValueError:
        return FALLBACK_COLOR


def base_group_color(base_group: str, base_groups: list[str]) -> str:
    """Color for a base group (EAA, Whey) using the GROUP_COLORS palette."""
    try:
        return GROUP_COLORS[sorted(base_groups).index(base_group) % len(GROUP_COLORS)]
    except ValueError:
        return FALLBACK_COLOR


# Private alias kept for modules that import the old underscore name.
_base_group_color = base_group_color


def taxon_color(taxon: str) -> str:
    """Return a consistent color for a taxon family name.

    Named taxa use the curated TAXA_COLORS palette.  Unknown taxa receive a
    deterministic color derived from the taxon name hash — the same taxon
    always renders with the same color across reports and server restarts,
    unlike the previous mutable-global accumulator.
    """
    if taxon in TAXA_COLORS:
        return TAXA_COLORS[taxon]
    return _EXTRA_TAXA_PALETTE[hash(taxon) % len(_EXTRA_TAXA_PALETTE)]


def hex_rgba(hex_color: str, alpha: float) -> str:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def _sorted_timepoints(rows: list[dict]) -> list[str]:
    """Thin wrapper — delegates to preprocessing.sorted_timepoints.

    Kept so callers that haven't updated their imports still work.
    Prefer importing sorted_timepoints from preprocessing directly.
    """
    from .preprocessing import sorted_timepoints
    return sorted_timepoints(rows)
