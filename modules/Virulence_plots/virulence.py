#!/usr/bin/env python3
"""
relative_abundance_plot.py
────────────────────────────────────────────────────────────────────────────────
Stacked Bar Plot — Relative Abundance per Sample
Designed to be wrapped in a Nextflow module.

USAGE:
    python3 relative_abundance_plot.py \
        --rel-abund  feature-table-relative-abundance.tsv \
        --metadata   sample-metadata.tsv \
        --out-plot   relative_abundance.png \
        [--taxonomy  taxonomy.tsv] \
        [--top-n     20] \
        [--min-pct   0.5]

NEXTFLOW PROCESS BLOCK:
    process REL_ABUND_PLOT {
        input:
        path rel_abund
        path metadata
        path taxonomy    // optional, pass '' to skip

        output:
        path "relative_abundance.png"

        script:
        \"\"\"
        python3 relative_abundance_plot.py \\
            --rel-abund ${rel_abund} \\
            --metadata  ${metadata} \\
            --taxonomy  ${taxonomy} \\
            --out-plot  relative_abundance.png
        \"\"\"
    }

INPUT FORMAT:
    rel-abund TSV : rows = Feature_ID (ASV/OTU hash), cols = sample IDs
    metadata  TSV : QIIME2-style; cols include 'sample-id', 'subject', etc.
                    The #q2:types directive row is skipped automatically.
    taxonomy  TSV : optional; cols 'Feature ID' and 'Taxon' (semicolon-delimited).
                    When provided, ASV hashes are replaced with readable genus/family
                    labels. Without it the top-N ASVs are labelled ASV_1..N.
────────────────────────────────────────────────────────────────────────────────
"""

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

warnings.filterwarnings('ignore')


# ══════════════════════════════════════════════════════════════════════════════
# ❶  VIRULENCE CONFIGURATION
#    ─────────────────────────
#    Any taxon whose label contains one of these keywords (case-insensitive)
#    is treated as virulent and coloured in the VIRULENT_REDS palette.
#    Add or remove entries freely — no other code changes needed.
# ══════════════════════════════════════════════════════════════════════════════
VIRULENT_KEYWORDS = [
    "Escherichia-Shigella",
    "Escherichia",
    "Shigella",
    "Enterococcus",
    "Proteobacteria",
    "Clostridioides",       # catches 'Clostridioides difficile'
    "Ruminococcus",
]

# Red shades for virulent taxa (darkest → brightest as more taxa are added)
VIRULENT_REDS = [
    "#FF0000", "#E8000D", "#CC0011", "#B50015",
    "#9E001A", "#870020", "#700025", "#59002B",
]


# ══════════════════════════════════════════════════════════════════════════════
# ❷  GENERAL CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

# Number of top ASVs to show individually; everything else → "Other"
TOP_N = 20

# ASVs with mean relative abundance below this % across all samples are always
# collapsed into "Other" regardless of TOP_N
MIN_MEAN_PCT = 0.5

# Metadata column used to order bars (set None to keep metadata file order)
ORDER_BY_COLUMN = "days-since-experiment-start"

# Blue/neutral palette for non-virulent taxa
BENIGN_PALETTE = [
    "#4E79A7", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1",
    "#FF9DA7", "#9C755F", "#1F77B4", "#2CA02C", "#9467BD",
    "#8C564B", "#E377C2", "#BCBD22", "#17BECF", "#AEC7E8",
    "#F28E2B", "#BAB0AC", "#6BAED6", "#74C476", "#9ECAE1",
]
OTHER_COLOUR = "#888888"

# Plot aesthetics
DARK_MODE = True
PLOT_DPI  = 160
BAR_WIDTH = 0.82
FIGSIZE_H = 8
MIN_FIG_W = 16


# ══════════════════════════════════════════════════════════════════════════════
# ❸  TAXONOMY HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def parse_label(taxon_str: str) -> str:
    """
    Extract the most specific readable label from a QIIME2 taxonomy string.
    Priority: genus > family > order > class > phylum.
    """
    rank_prefixes = [
        ('g__', ''), ('f__', ' (f)'), ('o__', ' (o)'),
        ('c__', ' (c)'), ('p__', ' (p)'),
    ]
    parts = [p.strip() for p in taxon_str.split(';')]
    for prefix, suffix in rank_prefixes:
        for part in reversed(parts):
            if part.startswith(prefix) and len(part) > len(prefix):
                return part[len(prefix):] + suffix
    return taxon_str.split(';')[-1].strip()


def load_taxonomy(path: str) -> dict:
    """Return {feature_id: readable_label} from a taxonomy TSV."""
    if not path or not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep='\t')
    df = df.rename(columns={df.columns[0]: 'Feature ID'})
    return {row['Feature ID']: parse_label(row['Taxon'])
            for _, row in df.iterrows()}


def is_virulent(label: str) -> bool:
    """Return True if the taxon label matches any virulence keyword."""
    lower = label.lower()
    return any(kw.lower() in lower for kw in VIRULENT_KEYWORDS)


# ══════════════════════════════════════════════════════════════════════════════
# ❹  DATA PREPARATION
# ══════════════════════════════════════════════════════════════════════════════

def load_metadata(path: str) -> pd.DataFrame:
    """Load QIIME2 metadata, dropping the #q2:types directive row."""
    df = pd.read_csv(path, sep='\t', dtype=str)
    df = df[df.iloc[:, 0] != '#q2:types'].reset_index(drop=True)
    df.columns = [c.strip() for c in df.columns]
    return df


def collapse_to_top_n(rel: pd.DataFrame, tax_map: dict,
                       top_n: int, min_pct: float) -> pd.DataFrame:
    """
    1. Rename Feature IDs using taxonomy labels (or ASV_N if no taxonomy).
    2. Aggregate rows with the same label.
    3. Keep top_n taxa by mean abundance; collapse the rest to 'Other'.
    Returns DataFrame: rows = taxon labels + 'Other', cols = sample IDs.
    """
    df = rel.copy()

    if tax_map:
        df.index = [tax_map.get(i, i) for i in df.index]
    else:
        rank = {fid: f"ASV_{i+1}"
                for i, fid in enumerate(
                    rel.mean(axis=1).sort_values(ascending=False).index)}
        df.index = [rank[i] for i in rel.index]

    df = df.groupby(level=0).sum()

    mean_abund = df.mean(axis=1)
    df = df[mean_abund >= min_pct]

    top_taxa = (mean_abund[mean_abund >= min_pct]
                .sort_values(ascending=False)
                .head(top_n).index.tolist())

    top_df    = df.loc[df.index.isin(top_taxa)]
    other_pct = 100.0 - top_df.sum(axis=0)
    other_row = pd.DataFrame([other_pct], index=['Other'])

    return pd.concat([top_df, other_row])


def order_samples(sample_ids: list, metadata: pd.DataFrame,
                  order_col: str) -> list:
    """Sort sample IDs by a numeric metadata column."""
    if not order_col or order_col not in metadata.columns:
        return sample_ids
    sub = metadata[metadata['sample-id'].isin(sample_ids)].copy()
    sub[order_col] = pd.to_numeric(sub[order_col], errors='coerce')
    sub = sub.sort_values(order_col)
    ordered = sub['sample-id'].tolist()
    ordered += [s for s in sample_ids if s not in ordered]
    return ordered


def build_colour_map(taxa: list) -> dict:
    """
    Assign colours to taxa:
      - Virulent taxa  → VIRULENT_REDS  (sequential red shades)
      - Benign taxa    → BENIGN_PALETTE (blue/neutral shades)
      - 'Other'        → OTHER_COLOUR   (grey)
    """
    colour_map = {}
    vir_idx = 0
    ben_idx = 0
    for taxon in taxa:
        if taxon == 'Other':
            colour_map[taxon] = OTHER_COLOUR
        elif is_virulent(taxon):
            colour_map[taxon] = VIRULENT_REDS[vir_idx % len(VIRULENT_REDS)]
            vir_idx += 1
        else:
            colour_map[taxon] = BENIGN_PALETTE[ben_idx % len(BENIGN_PALETTE)]
            ben_idx += 1
    colour_map['Other'] = OTHER_COLOUR
    return colour_map


# ══════════════════════════════════════════════════════════════════════════════
# ❺  PLOTTING
# ══════════════════════════════════════════════════════════════════════════════

BG    = '#0F1117'
PANEL = '#161B22'
FG    = '#CDD9E5'
GRID  = '#2A3040'
SPINE = '#444C56'


def _style(ax):
    ax.set_facecolor(PANEL if DARK_MODE else 'white')
    fg = FG if DARK_MODE else '#212529'
    ax.tick_params(colors=fg, labelsize=7.5)
    ax.xaxis.label.set_color(fg)
    ax.yaxis.label.set_color(fg)
    ax.title.set_color(fg)
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    for s in ['bottom', 'left']:
        ax.spines[s].set_color(GRID if DARK_MODE else '#DEE2E6')
    ax.yaxis.grid(True,
                  color=GRID if DARK_MODE else '#DEE2E6',
                  linewidth=0.6, linestyle='--')
    ax.set_axisbelow(True)


def draw_bars(ax, plot_df: pd.DataFrame, sample_ids: list,
              colour_map: dict, title: str):
    """
    Draw the stacked bars onto ax.
    Stack order: benign taxa (bottom) → virulent taxa → Other (top).
    This keeps virulent bars visually prominent at the top of each bar.
    """
    all_taxa = [t for t in plot_df.index if t != 'Other']
    benign_taxa   = [t for t in all_taxa if not is_virulent(t)]
    virulent_taxa = [t for t in all_taxa if is_virulent(t)]

    x      = np.arange(len(sample_ids))
    bottom = np.zeros(len(sample_ids))

    # ── Benign taxa (blues) — drawn first from the bottom ──
    for taxon in benign_taxa:
        vals = plot_df.loc[taxon, sample_ids].fillna(0).values
        ax.bar(x, vals, BAR_WIDTH, bottom=bottom,
               color=colour_map[taxon], zorder=3, linewidth=0)
        bottom += vals

    # ── Virulent taxa (reds) — drawn on top of benign ──
    for taxon in virulent_taxa:
        vals = plot_df.loc[taxon, sample_ids].fillna(0).values
        ax.bar(x, vals, BAR_WIDTH, bottom=bottom,
               color=colour_map[taxon], zorder=3, linewidth=0)
        bottom += vals

    # ── Other (grey) — drawn last at the very top ──
    vals = plot_df.loc['Other', sample_ids].fillna(0).values
    ax.bar(x, vals, BAR_WIDTH, bottom=bottom,
           color=OTHER_COLOUR, zorder=3, linewidth=0)

    ax.set_xticks(x)
    ax.set_xticklabels(sample_ids, rotation=55, ha='right', fontsize=7.5)
    ax.set_xlim(-0.6, len(sample_ids) - 0.4)
    ax.set_ylim(0, 105)
    ax.set_title(title, fontsize=11, fontweight='bold', pad=8)
    ax.set_ylabel('Relative Abundance (%)', fontsize=9)
    _style(ax)


def make_legend(fig, colour_map: dict, plot_df: pd.DataFrame):
    """
    Two-section legend: virulent taxa (red section) then benign taxa (blue section).
    """
    fg = FG if DARK_MODE else '#212529'

    all_taxa      = [t for t in plot_df.index if t != 'Other']
    virulent_taxa = [t for t in all_taxa if is_virulent(t)]
    benign_taxa   = [t for t in all_taxa if not is_virulent(t)]

    handles = []

    # Section header + virulent entries
    if virulent_taxa:
        handles.append(mpatches.Patch(
            color='none',
            label='── VIRULENT ──────────────'))
        for t in virulent_taxa:
            handles.append(mpatches.Patch(
                color=colour_map[t],
                label=t if len(t) <= 32 else t[:30] + '…'))

    # Section header + benign entries
    handles.append(mpatches.Patch(
        color='none',
        label='── OTHER TAXA ────────────'))
    for t in benign_taxa:
        handles.append(mpatches.Patch(
            color=colour_map[t],
            label=t if len(t) <= 32 else t[:30] + '…'))

    handles.append(mpatches.Patch(color=OTHER_COLOUR, label='Other (collapsed)'))

    fig.legend(
        handles=handles,
        loc='lower center',
        ncol=min(6, len(handles)),
        bbox_to_anchor=(0.5, -0.08),
        framealpha=0.15,
        facecolor=PANEL if DARK_MODE else 'white',
        edgecolor=SPINE,
        labelcolor=fg,
        fontsize=7.5,
    )


# ══════════════════════════════════════════════════════════════════════════════
# ❻  MAIN PLOT BUILDER
# ══════════════════════════════════════════════════════════════════════════════

def make_plot(plot_df: pd.DataFrame, metadata: pd.DataFrame,
              order_col: str, out_path: str):
    """Single-panel plot — all samples in one bar chart, ordered by order_col."""
    bg = BG if DARK_MODE else 'white'
    fg = FG if DARK_MODE else '#212529'

    all_samples = list(plot_df.columns)
    sample_ids  = order_samples(all_samples, metadata, order_col)

    # Build colour map once so legend and bars are consistent
    all_taxa   = list(plot_df.index)   # includes 'Other'
    colour_map = build_colour_map(all_taxa)

    # Count virulent taxa found
    virulent_found = [t for t in all_taxa
                      if t != 'Other' and is_virulent(t)]
    benign_found   = [t for t in all_taxa
                      if t != 'Other' and not is_virulent(t)]
    print(f"[INFO] Virulent taxa in top {TOP_N}: {virulent_found or 'none detected'}")
    print(f"[INFO] Benign  taxa in top {TOP_N}: {len(benign_found)}")

    # Figure sizing — scale width with number of samples
    fig_w = max(MIN_FIG_W, len(sample_ids) * 0.55 + 4)
    fig   = plt.figure(figsize=(fig_w, FIGSIZE_H), facecolor=bg)
    ax    = fig.add_subplot(111)

    draw_bars(ax, plot_df, sample_ids, colour_map,
              title='Relative Abundance per Sample')

    make_legend(fig, colour_map, plot_df)

    # Annotate virulence threshold line
    ax.axhline(80, color='#FF4444', linestyle=':', linewidth=1,
               alpha=0.5, zorder=4)
    ax.text(len(sample_ids) - 0.45, 81,
            '← high virulence zone', color='#FF6666',
            fontsize=7, ha='right', va='bottom', style='italic')

    fig.suptitle(
        f'Relative Abundance per Sample  '
        f'(top {TOP_N} taxa ≥{MIN_MEAN_PCT}% mean  |  '
        f'red = virulent, blue = other)',
        color=fg, fontsize=11, fontweight='bold', y=1.01
    )

    plt.savefig(out_path, dpi=PLOT_DPI, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"[INFO] Plot saved → {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
# ❼  CLI
# ══════════════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(
        description='Stacked bar plot of relative abundance per sample '
                    'with virulent taxa highlighted in red.')
    p.add_argument('--rel-abund', required=True,
                   help='Relative abundance TSV (rows=Feature_ID, cols=samples)')
    p.add_argument('--metadata',  required=True,
                   help='QIIME2 metadata TSV')
    p.add_argument('--taxonomy',  default=None,
                   help='Optional taxonomy TSV for readable labels')
    p.add_argument('--out-plot',  default='relative_abundance.png')
    p.add_argument('--top-n',     type=int,   default=TOP_N,
                   help=f'Top N taxa to display (default {TOP_N})')
    p.add_argument('--min-pct',   type=float, default=MIN_MEAN_PCT,
                   help=f'Min mean %% to show separately (default {MIN_MEAN_PCT})')
    p.add_argument('--order-by',  default=ORDER_BY_COLUMN,
                   help=f'Metadata column to order bars by (default: {ORDER_BY_COLUMN})')
    p.add_argument('--light',     action='store_true',
                   help='Use light background instead of dark')
    return p.parse_args()


def main():
    global TOP_N, MIN_MEAN_PCT, ORDER_BY_COLUMN, DARK_MODE

    args = parse_args()
    TOP_N          = args.top_n
    MIN_MEAN_PCT   = args.min_pct
    ORDER_BY_COLUMN = args.order_by
    DARK_MODE      = not args.light

    print("[INFO] Loading relative abundance table...")
    rel = pd.read_csv(args.rel_abund, sep='\t', index_col=0)
    print(f"       {rel.shape[0]} features × {rel.shape[1]} samples")

    print("[INFO] Loading metadata...")
    meta = load_metadata(args.metadata)
    print(f"       {len(meta)} samples in metadata")

    tax_map = {}
    if args.taxonomy:
        print("[INFO] Loading taxonomy...")
        tax_map = load_taxonomy(args.taxonomy)
        print(f"       {len(tax_map)} taxonomy entries loaded")
    else:
        print("[INFO] No taxonomy provided — ASVs labelled ASV_1, ASV_2, ...")
        print("       NOTE: virulence keyword matching requires taxonomy labels.")
        print("             Export taxonomy.qza from QIIME2 to enable red colouring.")

    print(f"[INFO] Collapsing to top {TOP_N} taxa (min mean {MIN_MEAN_PCT}%)...")
    plot_df = collapse_to_top_n(rel, tax_map, TOP_N, MIN_MEAN_PCT)
    print(f"       {len(plot_df) - 1} taxa shown + 'Other'")

    print("[INFO] Generating plot...")
    make_plot(plot_df, meta, ORDER_BY_COLUMN, args.out_plot)
    print("[INFO] Done.")


if __name__ == '__main__':
    main()
    