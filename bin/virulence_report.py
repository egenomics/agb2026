#!/usr/bin/env python3
import argparse
import os
import warnings
import sys

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

# Optional: set to an sra_id string to generate only that patient's plot.
# Leave as None to generate all 16 patient plots.
PATIENT_ID_FILTER = None    # e.g. "ERR1074192"


VIRULENT_GENERA = [
    "Escherichia-Shigella",
    "Enterococcus",
    "Clostridioides",
]



# Top N non-virulent genera shown individually; rest → "Other (benign)"
TOP_N_BENIGN = 10

# Metadata column linking ASV table columns to metadata rows
SAMPLE_ID_COL = "sra_id"

# Metadata column + value that identifies healthy samples
HEALTHY_COL   = "healthy"
HEALTHY_VALUE = "yes"


VIRULENT_RED  = "#FF0000"
BENIGN_BLUES  = [
    "#4E79A7", "#76B7B2", "#59A14F", "#EDC948",
    "#B07AA1", "#1F77B4", "#2CA02C", "#9467BD",
    "#17BECF", "#6BAED6",
]
OTHER_COLOUR  = "#AAAAAA"

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

    # ── Filter to samples present in metadata ─────────────────────────────────
    meta_samples = set(meta[SAMPLE_ID_COL].dropna())
    shared = [c for c in asv.columns if c in meta_samples]
    asv = asv[shared]
    print(f"[INFO] {len(shared)} samples matched between ASV table and metadata")

    # ── Relative abundance (computed after filtering so columns sum to 100%) ──
    col_totals = asv.sum(axis=0).replace(0, np.nan)
    rel = asv.div(col_totals, axis=1) * 100.0

    # ── Map ASV IDs → best available taxonomy label ───────────────────────────
    # Falls back through Genus → Family → Order → Class → Phylum → Unclassified
    # so that NaN genus values never cause groupby to silently drop rows.
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

    # ── Split healthy vs patient ───────────────────────────────────────────────
    healthy_ids = set(meta.loc[meta[HEALTHY_COL] == HEALTHY_VALUE, SAMPLE_ID_COL])
    patient_ids = set(meta.loc[meta[HEALTHY_COL] != HEALTHY_VALUE, SAMPLE_ID_COL])

    healthy_samples = [c for c in rel.columns if c in healthy_ids]
    patient_samples = [c for c in rel.columns if c in patient_ids]

    print(f"[INFO] Healthy samples  : {len(healthy_samples)}")
    print(f"[INFO] Patient samples  : {len(patient_samples)}")

    # Sanity check
    col_sums = rel.sum(axis=0).round(1)
    if not (col_sums == 100.0).all():
        print(f"[WARN] Column sum range: {col_sums.min()}–{col_sums.max()}%")
    else:
        print("[INFO] All columns sum to 100% ✓")

    return rel, patient_samples, healthy_samples, meta


def classify_taxa(rel):
    """
    Split relative abundance table into virulent and benign DataFrames.

    Returns:
      virulent_df  — rows are virulent genera
      benign_df    — rows are top N benign genera + 'Other (benign)'
      benign_colour — {genus: colour} dict (consistent across all patient plots)
    """
    all_taxa = rel.index.tolist()

    virulent_taxa = [
        t for t in all_taxa
        if any(v.lower() == t.lower() for v in VIRULENT_GENERA)
    ]
    benign_taxa = [t for t in all_taxa if t not in virulent_taxa]

    virulent_df = (rel.loc[virulent_taxa]
                   if virulent_taxa
                   else pd.DataFrame(index=[], columns=rel.columns))

    # Top N benign genera by mean abundance across ALL samples (consistent order)
    benign_rel  = rel.loc[benign_taxa]
    mean_abund  = benign_rel.mean(axis=1).sort_values(ascending=False)
    top_taxa    = mean_abund.head(TOP_N_BENIGN).index.tolist()
    top_benign  = benign_rel.loc[top_taxa]
    other       = benign_rel.loc[~benign_rel.index.isin(top_taxa)].sum(axis=0)
    other_row   = pd.DataFrame([other], index=["Other (benign)"])
    benign_df   = pd.concat([top_benign, other_row])

    # Colour map — fixed order so colours are the same in every patient plot
    benign_colour = {
        t: BENIGN_BLUES[i % len(BENIGN_BLUES)]
        for i, t in enumerate(top_taxa)
    }
    benign_colour["Other (benign)"] = OTHER_COLOUR

    print(f"[INFO] Virulent genera found: {virulent_taxa or 'none'}")

    return virulent_df, benign_df, benign_colour


def _style(ax):
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


def draw_bar(ax, virulent_df, benign_df, benign_colour,
             sample_ids, title, bar_width=0.45, show_ylabel=True):
    """
    Draw stacked bars for the given sample_ids onto ax.
    Stack order (bottom → top): benign genera → Other (benign) → virulent (red)
    """
    x      = np.arange(len(sample_ids))
    bottom = np.zeros(len(sample_ids))

    # Benign top genera
    for genus in [t for t in benign_df.index if t != "Other (benign)"]:
        vals = benign_df.loc[genus, sample_ids].fillna(0).values
        ax.bar(x, vals, bar_width, bottom=bottom,
               color=benign_colour.get(genus, OTHER_COLOUR),
               zorder=3, linewidth=0)
        bottom += vals

    # Other (benign)
    vals = benign_df.loc["Other (benign)", sample_ids].fillna(0).values
    ax.bar(x, vals, bar_width, bottom=bottom,
           color=OTHER_COLOUR, zorder=3, linewidth=0)
    bottom += vals

    # Virulent — single red segment grouping all virulent genera
    if not virulent_df.empty:
        vir_total = virulent_df.reindex(columns=sample_ids).fillna(0).sum(axis=0).values
        ax.bar(x, vir_total, bar_width, bottom=bottom,
               color=VIRULENT_RED, zorder=4, linewidth=0)

    ax.set_xticks(x)
    ax.set_xticklabels(sample_ids, rotation=30, ha='right', fontsize=8)
    ax.set_xlim(-0.6, len(sample_ids) - 0.4)
    ax.set_ylim(0, 108)
    ax.set_title(title, fontsize=10, fontweight='bold', pad=8)
    if show_ylabel:
        ax.set_ylabel("Relative Abundance (%)", fontsize=9)
    _style(ax)


def build_legend(fig, benign_colour, virulent_taxa_found):
    """Two-section legend: Virulent (red) then Benign (blues)."""
    handles = []

    handles.append(mpatches.Patch(color='none', label='━━ VIRULENT ━━━━━━━━'))
    vir_label = ", ".join(virulent_taxa_found) if virulent_taxa_found else "none detected"
    handles.append(mpatches.Patch(color=VIRULENT_RED,
                                  label=f"Virulent (grouped): {vir_label}"))

    handles.append(mpatches.Patch(color='none', label='━━ OTHER TAXA ━━━━━━'))
    for genus, colour in benign_colour.items():
        handles.append(mpatches.Patch(color=colour, label=genus))

    fig.legend(
        handles=handles,
        loc='lower center',
        ncol=4,
        bbox_to_anchor=(0.5, -0.12),
        framealpha=0.15,
        facecolor=PANEL,
        edgecolor=SPINE,
        labelcolor=FG,
        fontsize=7.5,
    )



def plot_patient(patient_id, virulent_df, benign_df, benign_colour,
                 healthy_samples, outdir):
    """
    Generate and save one plot for a single patient.
    Left panel  = patient sample bar
    Right panel = healthy cohort average bar
    """
    ctrl_label = f"Healthy avg\n(n={len(healthy_samples)})"

    # Build healthy control single-column DataFrames
    ctrl_virulent = pd.DataFrame(
        virulent_df[healthy_samples].mean(axis=1).values.reshape(-1, 1),
        index=virulent_df.index,
        columns=[ctrl_label]
    ) if not virulent_df.empty else pd.DataFrame(
        index=[], columns=[ctrl_label]
    )

    ctrl_benign = pd.DataFrame(
        benign_df[healthy_samples].mean(axis=1).values.reshape(-1, 1),
        index=benign_df.index,
        columns=[ctrl_label]
    )

    # ── Figure: 2 panels, patient wider than control ──────────────────────────
    fig, (ax_pat, ax_ctrl) = plt.subplots(
        1, 2,
        figsize=(9, 6),
        facecolor=BG,
        gridspec_kw={'width_ratios': [2, 1], 'wspace': 0.08}
    )

    # Patient bar
    draw_bar(ax_pat, virulent_df, benign_df, benign_colour,
             [patient_id],
             title=f"Patient: {patient_id}",
             bar_width=0.4,
             show_ylabel=True)

    # Healthy control bar
    draw_bar(ax_ctrl, ctrl_virulent, ctrl_benign, benign_colour,
             [ctrl_label],
             title="Healthy Control",
             bar_width=0.4,
             show_ylabel=False)

    # Shared vertical divider
    fig.add_artist(plt.Line2D(
        [ax_pat.get_position().x1 + 0.01,
         ax_pat.get_position().x1 + 0.01],
        [0.15, 0.92],
        transform=fig.transFigure,
        color=SPINE, linewidth=1.2, linestyle='--'
    ))

    # Legend + title
    virulent_taxa_found = list(virulent_df.index) if not virulent_df.empty else []
    build_legend(fig, benign_colour, virulent_taxa_found)

    fig.suptitle(
        f"Relative Abundance  |  Patient {patient_id}  |  Red = Virulent taxa",
        color=FG, fontsize=11, fontweight='bold', y=1.01
    )

    out_path = os.path.join(outdir, f"virulence_{patient_id}.png")
    plt.savefig(out_path, dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"[INFO] Saved → {out_path}")




def parse_args():
    p = argparse.ArgumentParser(
        description="Per-patient virulence relative abundance plots"
    )
    p.add_argument("--asv-table",   dest="asv_table",  default=None)
    p.add_argument("--taxonomy",    dest="taxonomy",   default=None)
    p.add_argument("--metadata",    dest="metadata",   default=None)
    p.add_argument("--outdir",      dest="outdir",     default=None)
    p.add_argument(
        "--patient-id", dest="patient_id", default=None,
        help="Optional: sra_id of a single patient to plot. "
             "If omitted, all non-healthy patients are plotted."
    )
    return p.parse_args()

## Mainnnnnnnnnn

if __name__ == '__main__':
    args = parse_args()

    # CLI args override section ❶ defaults
    if args.asv_table:  ASV_TABLE_PATH    = args.asv_table
    if args.taxonomy:   TAXONOMY_PATH     = args.taxonomy
    if args.metadata:   METADATA_PATH     = args.metadata
    if args.outdir:     OUTDIR            = args.outdir
    if args.patient_id: PATIENT_ID_FILTER = args.patient_id

    os.makedirs(OUTDIR, exist_ok=True)

    rel, patient_samples, healthy_samples, meta = load_data()
    virulent_df, benign_df, benign_colour = classify_taxa(rel)

    # ── Apply patient filter ──────────────────────────────────────────────────
    if PATIENT_ID_FILTER:
        if PATIENT_ID_FILTER not in patient_samples:
            raise ValueError(
                f"--patient-id '{PATIENT_ID_FILTER}' not found in non-healthy samples.\n"
                f"Available patient IDs: {patient_samples}"
            )
        targets = [PATIENT_ID_FILTER]
        print(f"[INFO] Generating plot for single patient: {PATIENT_ID_FILTER}")
    else:
        targets = patient_samples
        print(f"[INFO] Generating plots for all {len(targets)} patients...")

    # ── Generate one plot per patient ─────────────────────────────────────────
    for pid in targets:
        plot_patient(pid, virulent_df, benign_df, benign_colour,
                     healthy_samples, OUTDIR)

    print(f"\n[INFO] Done. {len(targets)} plot(s) saved to '{OUTDIR}/'")
