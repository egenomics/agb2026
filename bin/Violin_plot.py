#!/usr/bin/env python3

import argparse
import os
import warnings
import sys

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings('ignore')


ASV_TABLE_PATH  = sys.argv[1]
TAXONOMY_PATH   = sys.argv[2]     # --taxonomy
METADATA_PATH   = sys.argv[3]  # --metadata


VIRULENT_GENERA = [
    "Escherichia-Shigella",
    "Enterococcus",
    "Clostridioides",
]


SAMPLE_ID_COL = "sra_id"
HEALTHY_COL   = "healthy"
HEALTHY_VALUE = "yes"

# Colours
HEALTHY_COLOUR  = "#4A90D9"   # blue  — healthy group
PATIENT_COLOUR  = "#E07B5D"   # coral — non-healthy group
DOT_ALPHA       = 0.75        # opacity of individual patient dots
JITTER_STRENGTH = 0.06        # horizontal scatter of dots

# Dark theme
BG    = "#0F1117"
PANEL = "#161B22"
FG    = "#CDD9E5"
GRID  = "#2A3040"
SPINE = "#444C56"



def load_data():
    print("[INFO] Loading files...")
    asv  = pd.read_csv(ASV_TABLE_PATH, sep='\t', index_col=0)
    tax  = pd.read_csv(TAXONOMY_PATH,  sep='\t')
    meta = pd.read_csv(METADATA_PATH,  sep='\t')

    # Filter to samples in metadata, then compute relative abundance
    shared = [c for c in asv.columns if c in set(meta[SAMPLE_ID_COL].dropna())]
    asv = asv[shared]
    print(f"[INFO] {len(shared)} samples matched")

    col_totals = asv.sum(axis=0).replace(0, np.nan)
    rel = asv.div(col_totals, axis=1) * 100.0

    # Map ASV IDs to best taxonomy label
    def best_label(row):
        for col in ['Genus', 'Family', 'Order', 'Class', 'Phylum']:
            if col in row.index and pd.notna(row[col]) and str(row[col]).strip():
                suffix = '' if col == 'Genus' else f' ({col[0].lower()})'
                return str(row[col]) + suffix
        return 'Unclassified'

    tax['_label'] = tax.apply(best_label, axis=1)
    label_map = dict(zip(tax['ASV_ID'], tax['_label']))
    rel.index = [label_map.get(i, i) for i in rel.index]
    rel = rel.groupby(level=0).sum()

    # Split healthy / non-healthy
    healthy_ids = set(meta.loc[meta[HEALTHY_COL] == HEALTHY_VALUE, SAMPLE_ID_COL])
    patient_ids = set(meta.loc[meta[HEALTHY_COL] != HEALTHY_VALUE, SAMPLE_ID_COL])
    healthy_cols = [c for c in rel.columns if c in healthy_ids]
    patient_cols = [c for c in rel.columns if c in patient_ids]

    print(f"[INFO] Healthy n={len(healthy_cols)}, Non-healthy n={len(patient_cols)}")

    # Identify virulent genera present in this dataset
    virulent_found = [
        t for t in rel.index
        if any(v.lower() == t.lower() for v in VIRULENT_GENERA)
    ]
    print(f"[INFO] Virulent genera detected: {virulent_found}")

    if not virulent_found:
        raise ValueError(
            "No virulent genera found in the dataset. "
            "Check VIRULENT_GENERA list or taxonomy labels."
        )

    return rel, virulent_found, healthy_cols, patient_cols



def mann_whitney(a, b):
    """
    Mann-Whitney U test (non-parametric, appropriate for small/skewed samples).
    Returns (U statistic, p-value, significance label).
    """
    if len(a) < 3 or len(b) < 3:
        return None, None, "n too small"
    u_stat, p_val = stats.mannwhitneyu(a, b, alternative='two-sided')
    if p_val < 0.001:
        sig = "***"
    elif p_val < 0.01:
        sig = "**"
    elif p_val < 0.05:
        sig = "*"
    else:
        sig = "ns"
    return u_stat, p_val, sig



def _style(ax):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=FG, labelsize=9)
    ax.xaxis.label.set_color(FG)
    ax.yaxis.label.set_color(FG)
    ax.title.set_color(FG)
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    for s in ['bottom', 'left']:
        ax.spines[s].set_color(SPINE)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6, linestyle='--')
    ax.set_axisbelow(True)


def draw_violin(ax, healthy_vals, patient_vals, genus_name):

    np.random.seed(42)

    all_vals  = np.concatenate([healthy_vals, patient_vals])
    # Cap y-axis at 95th percentile; outliers still plotted as dots above
    y_cap     = np.percentile(all_vals[all_vals > 0], 95) if (all_vals > 0).any() else 1.0
    y_cap     = max(y_cap, 0.1)   # minimum visible range
    y_top     = y_cap * 1.45      # extra headroom for bracket + annotation

    groups = [
        (1, np.array(healthy_vals, dtype=float), HEALTHY_COLOUR),
        (2, np.array(patient_vals, dtype=float), PATIENT_COLOUR),
    ]

    for pos, vals, colour in groups:

        # Values clipped for violin shape (outliers excluded from KDE)
        vals_clipped = vals[vals <= y_cap]

        # ── Violin body ───────────────────────────────────────────────────────
        if len(np.unique(vals_clipped)) >= 3:
            parts = ax.violinplot(
                [vals_clipped], positions=[pos], widths=0.55,
                showmeans=False, showmedians=False, showextrema=False
            )
            for pc in parts['bodies']:
                pc.set_facecolor(colour)
                pc.set_edgecolor(colour)
                pc.set_alpha(0.45)
                pc.set_linewidth(1.2)
        else:
            # Too few unique values for KDE — draw a simple bar
            ax.bar(pos, np.median(vals), width=0.4,
                   color=colour, alpha=0.4, zorder=2)

        # ── IQR box (thick vertical line) ─────────────────────────────────────
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        q3_capped = min(q3, y_cap)
        ax.vlines(pos, q1, q3_capped, color=colour,
                  linewidth=7, alpha=0.85, zorder=3)

        # ── Median tick ───────────────────────────────────────────────────────
        med_plot = min(med, y_cap)
        ax.hlines(med_plot, pos - 0.12, pos + 0.12,
                  color='white', linewidth=2.2, zorder=4)

        # ── Jittered dots (all values, outliers capped for visibility) ─────────
        jitter = np.random.uniform(-JITTER_STRENGTH, JITTER_STRENGTH, len(vals))
        y_plot = np.clip(vals, 0, y_top * 0.95)   # cap dots within axes
        outlier_mask = vals > y_cap
        # Normal dots
        ax.scatter(
            np.full(len(vals), pos)[~outlier_mask] + jitter[~outlier_mask],
            y_plot[~outlier_mask],
            color='white', edgecolors=colour,
            s=28, linewidths=0.8, alpha=DOT_ALPHA, zorder=5
        )
        # Outlier dots — different marker (triangle) to signal clipping
        if outlier_mask.any():
            ax.scatter(
                np.full(outlier_mask.sum(), pos) + jitter[outlier_mask],
                np.full(outlier_mask.sum(), y_cap * 1.05),
                color=colour, marker='^', s=40,
                alpha=0.9, zorder=5,
                label=f"Outlier (>{y_cap:.1f}%)"
            )
            # Label the actual value next to each outlier triangle
            for xj, v in zip(jitter[outlier_mask], vals[outlier_mask]):
                ax.text(pos + xj + 0.07, y_cap * 1.06,
                        f"{v:.0f}%", color=colour,
                        fontsize=6.5, va='bottom')

        # ── Mean +/- SD — FIX: use axis transform so position is stable ────────
        ax.text(
            pos, 0.97,
            f"μ={vals.mean():.1f}%\nσ={vals.std():.1f}%",
            ha='center', va='top', color=FG, fontsize=7.5,
            transform=ax.get_xaxis_transform()   # x=data, y=axes fraction
        )

    # ── Significance bracket — FIX: use y_cap not raw max ────────────────────
    _, p_val, sig = mann_whitney(healthy_vals, patient_vals)
    y_sig     = y_cap * 1.15
    bracket_h = y_cap * 0.05

    ax.plot([1, 1, 2, 2],
            [y_sig - bracket_h, y_sig, y_sig, y_sig - bracket_h],
            color=FG, linewidth=0.9, zorder=6)
    p_label = f"p={p_val:.3f} {sig}" if p_val is not None else sig
    ax.text(1.5, y_sig + bracket_h * 0.3, p_label,
            ha='center', va='bottom', color=FG, fontsize=8.5, zorder=6)

    # ── Axes ──────────────────────────────────────────────────────────────────
    ax.set_xlim(0.4, 2.6)
    ax.set_ylim(-0.2, y_top)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(
        [f"Healthy\n(n={len(healthy_vals)})",
         f"Non-healthy\n(n={len(patient_vals)})"],
        fontsize=9
    )
    ax.set_title(f"{genus_name}", fontsize=11,
                 fontweight='bold', pad=10, color=FG, style='italic')
    ax.set_ylabel("Relative Abundance (%)", fontsize=9)

    # Clip line to show where y-axis was capped
    ax.axhline(y_cap, color=SPINE, linewidth=0.8,
               linestyle=':', alpha=0.7, zorder=1)
    ax.text(2.55, y_cap, f" 95th pct\n {y_cap:.1f}%",
            color=SPINE, fontsize=6, va='center')

    _style(ax)


def make_violin_plot(rel, virulent_found, healthy_cols, patient_cols, outdir):
    """
    Build the full figure — one subplot per virulent genus.
    Subplots are arranged in a single row.
    """
    n = len(virulent_found)
    fig_w = max(5 * n, 8)
    fig, axes = plt.subplots(
        1, n,
        figsize=(fig_w, 6),
        facecolor=BG,
        sharey=False     # each genus has its own y-scale — abundance ranges differ
    )
    if n == 1:
        axes = [axes]   # ensure iterable when only one genus

    for ax, genus in zip(axes, virulent_found):
        healthy_vals = rel.loc[genus, healthy_cols].fillna(0).values
        patient_vals = rel.loc[genus, patient_cols].fillna(0).values
        draw_violin(ax, healthy_vals, patient_vals, genus)

    # ── Shared legend ─────────────────────────────────────────────────────────
    legend_handles = [
        mpatches.Patch(color=HEALTHY_COLOUR, alpha=0.7, label="Healthy"),
        mpatches.Patch(color=PATIENT_COLOUR, alpha=0.7, label="Non-healthy"),
        plt.Line2D([0], [0], color='white', linewidth=2.5, label="Median"),
        plt.Line2D([0], [0], color=HEALTHY_COLOUR, linewidth=6,
                   alpha=0.85, label="IQR"),
        plt.scatter([], [], color='white', edgecolors='grey',
                    s=28, linewidths=0.8, label="Individual patient"),
    ]
    fig.legend(
        handles=legend_handles,
        loc='lower center',
        ncol=5,
        bbox_to_anchor=(0.5, -0.08),
        framealpha=0.15,
        facecolor=PANEL,
        edgecolor=SPINE,
        labelcolor=FG,
        fontsize=8.5,
    )

    fig.suptitle(
        "Virulent Genus Abundance  |  Healthy vs Non-healthy patients\n"
        "Mann-Whitney U test  (* p<0.05  ** p<0.01  *** p<0.001  ns = not significant)",
        color=FG, fontsize=11, fontweight='bold', y=1.03
    )

    plt.tight_layout()
    out_path = os.path.join(outdir, "virulence_violin.png")
    plt.savefig(out_path, dpi=160, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"[INFO] Saved → {out_path}")



def parse_args():
    p = argparse.ArgumentParser(
        description="Violin plot of virulent genus abundance: healthy vs non-healthy"
    )
    p.add_argument("--asv-table", dest="asv_table", default=None)
    p.add_argument("--taxonomy",  dest="taxonomy",  default=None)
    p.add_argument("--metadata",  dest="metadata",  default=None)
    p.add_argument(
        "--outdir", dest="outdir", default=None,
        help=(
            "Optional. Directory to save the output plot. "
            "Defaults to 'Virulence_analysis' in the current working directory. "
            "Created automatically if it does not exist."
        )
    )
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()

    if args.asv_table: ASV_TABLE_PATH = args.asv_table
    if args.taxonomy:  TAXONOMY_PATH  = args.taxonomy
    if args.metadata:  METADATA_PATH  = args.metadata

    OUTDIR = args.outdir if args.outdir else "Virulence_analysis"

    os.makedirs(OUTDIR, exist_ok=True)
    print(f"[INFO] Output directory: {os.path.abspath(OUTDIR)}")

    rel, virulent_found, healthy_cols, patient_cols = load_data()
    make_violin_plot(rel, virulent_found, healthy_cols, patient_cols, OUTDIR)
    print("[INFO] Done.")