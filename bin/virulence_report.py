#!/usr/bin/env python3

import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

warnings.filterwarnings('ignore')


ASV_TABLE_PATH  = sys.argv[1]
TAXONOMY_PATH   = sys.argv[2]     # --taxonomy
METADATA_PATH   = sys.argv[3]  # --metadata
OUTPUT_PATH     = "virulence_abundance_plot.png"  # --output


# Exact genus-level labels as they appear in Group B's taxonomy output.
# "Escherichia-Shigella" is a single merged label (SILVA/QIIME2 convention).
VIRULENT_GENERA = [
    "Escherichia-Shigella",
    "Enterococcus",
    "Clostridioides",
]


TOP_N_BENIGN = 10

# Matches the actual column header in sample-metadata.tsv
METADATA_SAMPLE_ID_COL = "sample-id"

# 'healthy' column uses lowercase "yes" / "no"
HEALTHY_COL   = "healthy"
HEALTHY_VALUE = "yes"


VIRULENT_REDS = [
    "#FF0000", "#D10000", "#A80000",
    "#800000", "#FF4444", "#FF7777",
]

BENIGN_BLUES = [
    "#4E79A7", "#76B7B2", "#59A14F", "#EDC948",
    "#B07AA1", "#1F77B4", "#2CA02C", "#9467BD",
    "#17BECF", "#6BAED6",
]

OTHER_COLOUR   = "#AAAAAA"
VIRULENT_GROUP = "#FF0000"

BG    = "#0F1117"
PANEL = "#161B22"
FG    = "#CDD9E5"
GRID  = "#2A3040"
SPINE = "#444C56"


def load_and_process():
    print("[INFO] Loading files...")
    asv  = pd.read_csv(ASV_TABLE_PATH,  sep='\t', index_col=0)
    tax  = pd.read_csv(TAXONOMY_PATH,   sep='\t')
    meta = pd.read_csv(METADATA_PATH,   sep='\t')

    meta_samples = set(meta[METADATA_SAMPLE_ID_COL].dropna())
    shared = [c for c in asv.columns if c in meta_samples]
    asv = asv[shared]
    print(f"[INFO] {len(shared)} samples matched between ASV table and metadata")

    print("[INFO] Computing relative abundance...")
    col_totals = asv.sum(axis=0).replace(0, np.nan)
    rel = asv.div(col_totals, axis=1) * 100.0

    def best_label(row):
        for col in ['Genus', 'Family', 'Order', 'Class', 'Phylum']:
            if col in row and pd.notna(row[col]) and str(row[col]).strip():
                suffix = '' if col == 'Genus' else f' ({col[0].lower()})'
                return str(row[col]) + suffix
        return 'Unclassified'

    tax['_label'] = tax.apply(best_label, axis=1)
    genus_map = dict(zip(tax['ASV_ID'], tax['_label']))
    rel.index = [genus_map.get(i, i) for i in rel.index]
    rel = rel.groupby(level=0).sum()

    healthy_ids = set(
        meta.loc[meta[HEALTHY_COL] == HEALTHY_VALUE, METADATA_SAMPLE_ID_COL]
    )
    patient_ids = set(
        meta.loc[meta[HEALTHY_COL] != HEALTHY_VALUE, METADATA_SAMPLE_ID_COL]
    )

    healthy_samples = [c for c in rel.columns if c in healthy_ids]
    patient_samples = [c for c in rel.columns if c in patient_ids]

    print(f"[INFO] Healthy samples: {len(healthy_samples)}")
    print(f"[INFO] Patient (non-healthy) samples: {len(patient_samples)}")

    return rel, patient_samples, healthy_samples


def classify_taxa(rel):
    all_genera = rel.index.tolist()

    virulent_genera = [
        g for g in all_genera
        if any(v.lower() == g.lower() for v in VIRULENT_GENERA)
    ]
    benign_genera = [g for g in all_genera if g not in virulent_genera]

    virulent_df = rel.loc[virulent_genera] if virulent_genera else pd.DataFrame(columns=rel.columns)

    benign_rel   = rel.loc[benign_genera]
    mean_abund   = benign_rel.mean(axis=1).sort_values(ascending=False)
    top_genera   = mean_abund.head(TOP_N_BENIGN).index.tolist()
    top_benign   = benign_rel.loc[top_genera]

    other_benign = benign_rel.loc[~benign_rel.index.isin(top_genera)].sum(axis=0)
    other_row    = pd.DataFrame([other_benign], index=["Other (benign)"])
    benign_df    = pd.concat([top_benign, other_row])

    total = virulent_df.sum(axis=0) if not virulent_df.empty else 0
    total = total + benign_df.sum(axis=0)
    if not (total.round(1) == 100.0).all():
        print(f"[WARN] Some columns don't sum to 100% — min={total.min():.1f}% max={total.max():.1f}%")
    else:
        print("[INFO] All sample columns sum to 100% ✓")

    print(f"[INFO] Virulent genera detected: {virulent_genera or 'none'}")
    print(f"[INFO] Top {TOP_N_BENIGN} named benign genera + 'Other (benign)' shown")

    return virulent_df, benign_df


def _style_ax(ax):
    ax.set_facecolor(PANEL)
    ax.tick_params(colors=FG)
    ax.xaxis.label.set_color(FG)
    ax.yaxis.label.set_color(FG)
    ax.title.set_color(FG)
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    for s in ['bottom', 'left']:
        ax.spines[s].set_color(SPINE)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6, linestyle='--')
    ax.set_axisbelow(True)


def draw_stacked_bars(ax, virulent_df, benign_df, sample_ids,
                      title, show_ylabel=True, is_control=False):
    x      = np.arange(len(sample_ids))
    bottom = np.zeros(len(sample_ids))
    width  = 0.65 if not is_control else 0.35

    benign_genera = [g for g in benign_df.index if g != "Other (benign)"]
    benign_colour = {g: BENIGN_BLUES[i % len(BENIGN_BLUES)]
                     for i, g in enumerate(benign_genera)}

    for genus in benign_genera:
        vals = benign_df.loc[genus, sample_ids].fillna(0).values
        ax.bar(x, vals, width, bottom=bottom,
               color=benign_colour[genus], zorder=3, linewidth=0)
        bottom += vals

    if "Other (benign)" in benign_df.index:
        vals = benign_df.loc["Other (benign)", sample_ids].fillna(0).values
        ax.bar(x, vals, width, bottom=bottom,
               color=OTHER_COLOUR, zorder=3, linewidth=0)
        bottom += vals

    if not virulent_df.empty:
        virulent_total = virulent_df[sample_ids].fillna(0).sum(axis=0).values
        ax.bar(x, virulent_total, width, bottom=bottom,
               color=VIRULENT_GROUP, zorder=4, linewidth=0,
               label="Virulent (grouped)")
        bottom += virulent_total

    ax.set_xticks(x)
    ax.set_xticklabels(
        sample_ids, rotation=60, ha='right',
        fontsize=6.5 if not is_control else 9
    )
    ax.set_xlim(-0.6, len(sample_ids) - 0.4)
    ax.set_ylim(0, 108)
    ax.set_title(title, fontsize=10, fontweight='bold', pad=8, color=FG)
    if show_ylabel:
        ax.set_ylabel("Relative Abundance (%)", fontsize=9, color=FG)
    _style_ax(ax)

    return benign_colour


def build_legend(fig, benign_colour, virulent_genera_found):
    handles = []

    handles.append(mpatches.Patch(color='none', label='━━ VIRULENT TAXA ━━'))
    label = (", ".join(virulent_genera_found)
             if virulent_genera_found else "none detected")
    handles.append(mpatches.Patch(color=VIRULENT_GROUP,
                                  label=f"Virulent group: {label}"))

    handles.append(mpatches.Patch(color='none', label='━━ OTHER TAXA ━━━━'))
    for genus, col in benign_colour.items():
        handles.append(mpatches.Patch(color=col, label=genus))
    handles.append(mpatches.Patch(color=OTHER_COLOUR, label="Other (benign)"))

    fig.legend(
        handles=handles,
        loc='lower center',
        ncol=4,
        bbox_to_anchor=(0.5, -0.1),
        framealpha=0.15,
        facecolor=PANEL,
        edgecolor=SPINE,
        labelcolor=FG,
        fontsize=7.5,
    )


def make_plot(rel, patient_samples, healthy_samples):
    virulent_df, benign_df = classify_taxa(rel)

    healthy_virulent_mean = (
        virulent_df[healthy_samples].mean(axis=1)
        if not virulent_df.empty and healthy_samples
        else pd.Series(dtype=float)
    )
    healthy_benign_mean = benign_df[healthy_samples].mean(axis=1)

    ctrl_label = f"Healthy\nControl\n(n={len(healthy_samples)})"
    ctrl_virulent = pd.DataFrame(
        healthy_virulent_mean.values.reshape(-1, 1),
        index=healthy_virulent_mean.index,
        columns=[ctrl_label]
    ) if not healthy_virulent_mean.empty else pd.DataFrame(columns=[ctrl_label])

    ctrl_benign = pd.DataFrame(
        healthy_benign_mean.values.reshape(-1, 1),
        index=healthy_benign_mean.index,
        columns=[ctrl_label]
    )

    n_patients = len(patient_samples)
    fig_w = max(16, n_patients * 0.45 + 5)
    fig   = plt.figure(figsize=(fig_w, 8), facecolor=BG)

    gs = fig.add_gridspec(
        1, 2,
        width_ratios=[n_patients, 2],
        wspace=0.06
    )
    ax_patients = fig.add_subplot(gs[0])
    ax_control  = fig.add_subplot(gs[1])

    benign_colour = draw_stacked_bars(
        ax_patients, virulent_df, benign_df,
        patient_samples,
        title=f"Patient Samples  (n={len(patient_samples)})",
        show_ylabel=True,
        is_control=False
    )

    draw_stacked_bars(
        ax_control, ctrl_virulent, ctrl_benign,
        [ctrl_label],
        title="Healthy\nControl Avg",
        show_ylabel=False,
        is_control=True
    )

    fig.add_artist(
        plt.Line2D(
            [ax_patients.get_position().x1 + 0.005,
             ax_patients.get_position().x1 + 0.005],
            [0.1, 0.92],
            transform=fig.transFigure,
            color=SPINE, linewidth=1.2, linestyle='--'
        )
    )

    virulent_found = [
        g for g in rel.index
        if any(v.lower() == g.lower() for v in VIRULENT_GENERA)
    ]
    build_legend(fig, benign_colour, virulent_found)

    fig.suptitle(
        "Relative Abundance per Sample  |  Red = Virulent taxa (grouped)",
        color=FG, fontsize=12, fontweight='bold', y=1.01
    )

    plt.savefig(OUTPUT_PATH, dpi=160, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"[INFO] Plot saved → {OUTPUT_PATH}")


def parse_args():
    import argparse
    p = argparse.ArgumentParser(
        description="Virulence relative abundance stacked bar plot"
    )
    p.add_argument("--asv-table", dest="asv_table", default=None)
    p.add_argument("--taxonomy",  dest="taxonomy",  default=None)
    p.add_argument("--metadata",  dest="metadata",  default=None)
    p.add_argument("--output",    dest="output",    default=None)
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()

    if args.asv_table: ASV_TABLE_PATH = args.asv_table
    if args.taxonomy:  TAXONOMY_PATH  = args.taxonomy
    if args.metadata:  METADATA_PATH  = args.metadata
    if args.output:    OUTPUT_PATH    = args.output

    rel, patient_samples, healthy_samples = load_and_process()
    make_plot(rel, patient_samples, healthy_samples)
    print("[INFO] Done.")