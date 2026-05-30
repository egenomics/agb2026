#!/usr/bin/env python3
"""
Picks a rarefaction depth automatically from QIIME2 rarefaction curves.

Steps:
  1. Unzip the .qzv and grab the per-sample curve CSV/TSVs
  2. Find where each sample's curve flattens
  3. Pick a global depth that covers enough samples
  4. Save plots + output files
  
Methods:
  coverage
      Chooses the depth at which X% of samples have reached a plateau. Recommended.
  percentile
      Uses the Nth percentile of the per-sample plateau depths.
  knee
      Detects the knee point from each individual sample curve (not from the overall distribution).
      The global threshold is then set to the median plateau depth across samples that pass QC.
      
Usage:
    python select_rarefaction_depth.py
        --curves      rarefaction_curves.qzv
        --shannon     diversity_table/shannon/alpha-diversity.tsv
        --observed    diversity_table/observed/alpha-diversity.tsv
        --faith       diversity_table/faith/alpha-diversity.tsv
        --simpson     diversity_table/simpson/alpha-diversity.tsv
        --output      rarefaction_output/
        --method      coverage
        --coverage-pct 90
        --metric      observed_features
        --dropout-max 0.10
"""

import argparse
import json
import logging
import warnings
import zipfile
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

try:
    import matplotlib
    matplotlib.use("Agg")  # no display needed
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    from matplotlib.lines import Line2D
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False
    log.warning("matplotlib not found, skipping plots")


# step 1: get the curves out of the .qzv (it's just a zip file)

def extract_curves_from_qzv(qzv_path, metric):
    log.info("Extracting curves from %s (metric: %s)", qzv_path, metric)

    with zipfile.ZipFile(qzv_path, "r") as zf:
        candidates = [n for n in zf.namelist()
                      if n.endswith(f"{metric}.csv") and "/data/" in n]
        if not candidates:
            candidates = [n for n in zf.namelist()
                          if metric in n and n.endswith(".csv")]
        if not candidates:
            available = [n for n in zf.namelist() if n.endswith(".csv")]
            raise FileNotFoundError(
                f"Can't find '{metric}.csv' in {qzv_path}.\n"
                f"CSVs in archive: {available}")
                
        with zf.open(candidates[0]) as f:
            raw = pd.read_csv(f, index_col=0)

    # columns are like "depth-500_iter-1" parse depth and average over iterations
    depth_cols = {}
    for col in raw.columns:
        col_clean = col.lower().replace("depth-", "").replace("depth_", "")
        for sep in ["_iter-", "_iter_"]:
            if sep in col_clean:
                try:
                    depth = int(col_clean.split(sep)[0])
                    depth_cols.setdefault(depth, []).append(col)
                except ValueError:
                    pass
                break

    if not depth_cols:
        raise ValueError(f"Couldn't parse depth columns. First few: {list(raw.columns[:5])}")

    curves = {d: raw[cols].mean(axis=1) for d, cols in sorted(depth_cols.items())}
    curve_df = pd.DataFrame(curves)
    curve_df.index.name = "sample-id"

    log.info("%d samples x %d depths (%d to %d reads)", *curve_df.shape, curve_df.columns.min(), curve_df.columns.max())
    return curve_df


# step 2: find where each sample's curve levels off

def find_knee_on_curve(depths, diversity):
    # knee method: find the point furthest from the line connecting start and end
    valid = ~np.isnan(diversity)
    x = depths[valid].astype(float)
    y = diversity[valid].astype(float)

    if len(x) < 3:
        return float(x[-1]) if len(x) > 0 else np.nan

    # normalize to [0,1] so x and y are on the same scale
    x_n = (x - x[0]) / (x[-1] - x[0] + 1e-12)
    y_n = (y - y.min()) / (y.max() - y.min() + 1e-12)

    dx, dy = x_n[-1] - x_n[0], y_n[-1] - y_n[0]
    dists = np.abs(dy * x_n - dx * y_n + x_n[-1] * y_n[0] - y_n[-1] * x_n[0]) / (np.sqrt(dx**2 + dy**2) + 1e-12)

    return float(x[int(np.argmax(dists))])


def per_sample_plateau_depths(curve_df):
    depths = curve_df.columns.values.astype(float)
    result = {}
    for sample_id, row in curve_df.iterrows():
        y = row.values.astype(float)
        if np.sum(~np.isnan(y)) < 3:
            result[sample_id] = np.nan
        else:
            result[sample_id] = find_knee_on_curve(depths, y)

    series = pd.Series(result, name="plateau_depth")
    log.info("Plateaus: min=%.0f  median=%.0f  max=%.0f  NaN=%d",
             series.dropna().min(), series.dropna().median(), series.dropna().max(), series.isna().sum())
    return series


# step 3: pick one global threshold from the per-sample plateaus

def select_global_threshold(plateau_depths, method, percentile, coverage_pct, dropout_max):
    
    # a sample "passes" if its curve already flattened at or before the threshold
    valid = plateau_depths.dropna()
    n_total = len(plateau_depths)
    n_nan = int(plateau_depths.isna().sum())

    log.info("Step 3/5 Selecting threshold (method=%s)...", method)

    if method == "coverage":
        # depth where coverage_pct% of samples have already plateaued
        threshold = int(np.round(np.percentile(valid.values, coverage_pct)))
    elif method == "percentile":
        threshold = int(np.round(np.percentile(valid.values, percentile)))
    elif method == "knee":
        threshold = int(np.round(np.median(valid.values)))
    else:
        raise ValueError(f"Unknown method: {method!r}")

    log.info("threshold = %d reads", threshold)

    pass_mask = plateau_depths <= threshold
    pass_mask[plateau_depths.isna()] = False  # too shallow = always excluded

    n_pass = int(pass_mask.sum())
    n_fail = n_total - n_pass
    dropout_frac = n_fail / n_total

    log.info("Retained %d / %d (dropped %.1f%%)", n_pass, n_total, dropout_frac * 100)

    if n_pass == 0:
        warnings.warn(
            f"No samples retained at threshold={threshold}. "
            "Try --method coverage --coverage-pct 90 or increase --percentile."
        )
    if dropout_frac > dropout_max:
        warnings.warn(f"Dropout {dropout_frac:.1%} exceeds --dropout-max {dropout_max:.1%}.")

    return {
        "threshold":      threshold,
        "method":         method,
        "pass_mask":      pass_mask,
        "plateau_depths": plateau_depths,
        "n_pass":         n_pass,
        "n_fail":         n_fail,
        "n_nan":          n_nan,
        "dropout_frac":   dropout_frac,
    }


# step 4: load the alpha diversity tables from module 1

def load_alpha_tsv(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df.iloc[:, 0].rename(df.columns[0])


def load_module1_metrics(shannon, observed, faith, simpson):
    series = {}
    for name, path in [("shannon", shannon), ("observed_features", observed), ("faith_pd", faith), ("simpson", simpson)]:
        if path and Path(path).exists():
            try:
                series[name] = load_alpha_tsv(path)
                log.info("Loaded %s (%d samples)", name, len(series[name]))
            except Exception as e:
                log.warning("Could not load %s: %s", name, e)
    
    return pd.DataFrame(series) if series else pd.DataFrame()


# step 5: plots

def plot_rarefaction_curves(curve_df, pass_mask, plateau_depths, threshold, out_dir):
    fig, ax = plt.subplots(figsize=(13, 6))
    depths = curve_df.columns.values.astype(float)

    for sample_id, row in curve_df.iterrows():
        retained = bool(pass_mask.get(sample_id, False))
        colour = "#2980B9" if retained else "#E74C3C"
        alpha  = 0.65      if retained else 0.20
        ax.plot(depths, row.values, color=colour, alpha=alpha, linewidth=0.8)

        # dot at the detected plateau point
        pd_val = plateau_depths.get(sample_id, np.nan)
        if not np.isnan(pd_val):
            closest = depths[np.argmin(np.abs(depths - pd_val))]
            ax.scatter(closest, row[closest], color=colour, s=20, alpha=0.6, zorder=3)

    ax.axvline(threshold, color="#2C3E50", linestyle="--", linewidth=2)
    
    n_ret = int(pass_mask.sum())
    n_exc = len(pass_mask) - n_ret
    handles = [
        Line2D([0], [0], color="#2980B9", lw=2, label=f"Retained (n={n_ret})"),
        Line2D([0], [0], color="#E74C3C", lw=2, label=f"Excluded (n={n_exc})"),
        Line2D([0], [0], color="#2C3E50", lw=2, linestyle="--",
               label=f"Threshold = {threshold:,} reads"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=9)
    ax.set_xlabel("Sequencing depth (reads)", fontsize=11)
    ax.set_ylabel("Alpha diversity", fontsize=11)
    ax.set_title("Rarefaction curves (dot = per-sample plateau, dashed = global threshold)", fontsize=11)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    fig.tight_layout()
    fig.savefig(out_dir / "rarefaction_curves.png", dpi=150)
    plt.close(fig)
    log.info("Saved rarefaction_curves.png")


def plot_plateau_distribution(result, out_dir):
    plateau_depths = result["plateau_depths"].dropna()
    threshold = result["threshold"]
    pass_mask = result["pass_mask"]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # left: one bar per sample sorted by plateau depth
    ax = axes[0]
    sorted_pd = np.sort(plateau_depths.values)
    colours = ["#2980B9" if v <= threshold else "#E74C3C" for v in sorted_pd]
    ax.bar(range(len(sorted_pd)), sorted_pd, color=colours, alpha=0.85, width=1.0)
    ax.axhline(threshold, color="#2C3E50", linestyle="--", linewidth=2,
               label=f"Threshold = {threshold:,} reads")
    ax.set_xlabel("Samples (sorted)")
    ax.set_ylabel("Plateau depth (reads)")
    ax.set_title("Per-sample plateau depths (blue = retained)")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.legend(frameon=False)

    # right: histogram
    ax2 = axes[1]
    pm_aligned = pass_mask.reindex(plateau_depths.index).fillna(False)
    retained_pd = plateau_depths[pm_aligned]
    excluded_pd = plateau_depths[~pm_aligned]
    bins = np.linspace(plateau_depths.min(), plateau_depths.max(), 25)
    
    if len(retained_pd):
        ax2.hist(retained_pd.values, bins=bins, color="#2980B9", alpha=0.75,
                 label=f"Retained (n={len(retained_pd)})")
                 
    if len(excluded_pd):
        ax2.hist(excluded_pd.values, bins=bins, color="#E74C3C", alpha=0.75,
                 label=f"Excluded (n={len(excluded_pd)})")
    ax2.axvline(threshold, color="#2C3E50", linestyle="--", linewidth=2,
                label=f"Threshold = {threshold:,} reads")
    ax2.set_xlabel("Plateau depth (reads)")
    ax2.set_ylabel("Number of samples")
    ax2.set_title("Distribution of plateau depths")
    ax2.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax2.legend(frameon=False)

    fig.suptitle(
        f"{result['method']} method -> threshold = {threshold:,} reads  "
        f"(retained: {result['n_pass']}, excluded: {result['n_fail']})",
        fontsize=11
    )
    fig.tight_layout()
    fig.savefig(out_dir / "plateau_depth_distribution.png", dpi=150)
    plt.close(fig)
    log.info("Saved plateau_depth_distribution.png")


def plot_module1_diversity(diversity_df, pass_mask, threshold, out_dir):
    metrics = [c for c in ["observed_features", "shannon", "faith_pd", "simpson"]
               if c in diversity_df.columns]
    if not metrics:
        return

    fig, axes = plt.subplots(1, len(metrics), figsize=(4.5 * len(metrics), 5), squeeze=False)
    for i, metric in enumerate(metrics):
        ax = axes[0][i]
        vals = diversity_df[metric].dropna()
        common = vals.index.intersection(pass_mask.index)
        ret = vals.loc[common][pass_mask.loc[common]]
        exc = vals.loc[common][~pass_mask.loc[common]]

        parts, labels, colours = [], [], []
        if len(ret):
            parts.append(ret.values)
            labels.append(f"Retained\n(n={len(ret)})")
            colours.append("#2980B9")
        if len(exc):
            parts.append(exc.values)
            labels.append(f"Excluded\n(n={len(exc)})")
            colours.append("#E74C3C")

        if parts:
            vp = ax.violinplot(parts, positions=range(len(parts)), showmedians=True, showextrema=True)
            for pc, col in zip(vp["bodies"], colours):
                pc.set_facecolor(col)
                pc.set_alpha(0.7)
            for part in ["cmedians", "cmaxes", "cmins", "cbars"]:
                if part in vp:
                    vp[part].set_color("#2C3E50")

        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel(metric.replace("_", " ").title())
        ax.set_title(metric.replace("_", " ").title())

    fig.suptitle(f"Module 1 alpha diversity split by retention (threshold = {threshold:,} reads)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_dir / "module1_diversity_by_retention.png", dpi=150)
    plt.close(fig)
    log.info("Saved module1_diversity_by_retention.png")


# main

def main():
    parser = argparse.ArgumentParser(
        description="Pick rarefaction depth from QIIME2 curves.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        
    parser.add_argument("--curves", required=True)
    parser.add_argument("--shannon", default=None)
    parser.add_argument("--observed", default=None)
    parser.add_argument("--faith", default=None)
    parser.add_argument("--simpson", default=None)
    parser.add_argument("--output", default="rarefaction_output")
    parser.add_argument("--method", default="coverage", choices=["coverage", "percentile", "knee"])
    parser.add_argument("--metric", default="observed_features")
    parser.add_argument("--coverage-pct", type=float, default=90, dest="coverage_pct", help="keep the %% of samples with lowest plateau depths")
    parser.add_argument("--percentile", type=float, default=75)
    parser.add_argument("--dropout-max", type=float, default=0.10, dest="dropout_max")
    parser.add_argument("--sampling-depth", type=int,   default=1103, dest="sampling_depth", help="module 1 depth, only used for plot labels")
    
    args = parser.parse_args()

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("Step 1/5 Extracting rarefaction curves...")
    curve_df = extract_curves_from_qzv(args.curves, args.metric)

    log.info("Step 2/5 Finding plateau depth for each sample...")
    
    plateau_depths = per_sample_plateau_depths(curve_df)

    result = select_global_threshold(plateau_depths, args.method, args.percentile, args.coverage_pct, args.dropout_max)

    log.info("Step 4/5 Loading alpha diversity tables...")
    diversity_df = load_module1_metrics(args.shannon, args.observed, args.faith, args.simpson)

    log.info("Step 5/5 Saving outputs...")
    
    if HAS_PLOT:
        plot_rarefaction_curves(curve_df, result["pass_mask"], plateau_depths, result["threshold"], out_dir)
        plot_plateau_distribution(result, out_dir)
        if not diversity_df.empty:
            plot_module1_diversity(diversity_df, result["pass_mask"], result["threshold"], out_dir)

    (out_dir / "rarefaction_threshold.txt").write_text(str(result["threshold"]) + "\n")

    keep = result["pass_mask"].index[result["pass_mask"]].tolist()
    pd.DataFrame({"sample-id": keep}).to_csv(out_dir / "samples_to_keep.tsv", sep="\t", index=False)

    qc_df = pd.DataFrame({
        "sample_id":     plateau_depths.index,
        "plateau_depth": plateau_depths.values,
        "passes":        result["pass_mask"].values,
    })
    
    if not diversity_df.empty:
        for col in ["observed_features", "shannon", "faith_pd", "simpson"]:
            if col in diversity_df.columns:
                qc_df = qc_df.merge(diversity_df[[col]].rename_axis("sample_id").reset_index(),on="sample_id", how="left")
    qc_df.to_csv(out_dir / "sample_qc.tsv", sep="\t", index=False)

    summary = {
        "date": datetime.now().isoformat(),
        "method": args.method,
        "curve_metric": args.metric,
        "threshold": result["threshold"],
        "n_total": len(plateau_depths),
        "n_retained": result["n_pass"],
        "n_dropped": result["n_fail"],
        "n_too_shallow": result["n_nan"],
        "dropout_frac": round(result["dropout_frac"], 4),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    report = [
        "=== Rarefaction threshold selection ===",
        f"Date      : {summary['date']}",
        f"Method    : {args.method}",
        f"Metric    : {args.metric}",
        f"Threshold : {result['threshold']:,} reads",
        f"Total     : {len(plateau_depths)} samples",
        f"Retained  : {result['n_pass']}",
        f"Excluded  : {result['n_fail']} ({result['dropout_frac']:.1%})",
        f"  too shallow: {result['n_nan']}",
        "",
        "Output files:",
        "  rarefaction_threshold.txt          -- use this in qiime feature-table rarefy",
        "  samples_to_keep.tsv                -- sample list for qiime filter-samples",
        "  sample_qc.tsv                      -- per-sample plateau depth + pass/fail",
        "  summary.json",
        "  rarefaction_curves.png",
        "  plateau_depth_distribution.png",
        "  module1_diversity_by_retention.png",
    ]
    (out_dir / "report.txt").write_text("\n".join(report) + "\n")
    print("\n".join(report))


    # merge the 3 individual PNGs into one stacked image
    individual_pngs = [
        out_dir / "rarefaction_curves.png",
        out_dir / "plateau_depth_distribution.png",
        out_dir / "module1_diversity_by_retention.png",
    ]
    existing = [p for p in individual_pngs if p.exists()]
    if existing:
        imgs = [plt.imread(str(p)) for p in existing]
        fig, axes = plt.subplots(len(imgs), 1,
                                 figsize=(13, 6 * len(imgs)),
                                 constrained_layout=True)
        if len(imgs) == 1:
            axes = [axes]
        titles = [p.stem.replace("_", " ").title() for p in existing]
        for ax, img, title in zip(axes, imgs, titles):
            ax.imshow(img)
            ax.axis("off")
            ax.set_title(title, fontsize=13, pad=8)
        fig.savefig(out_dir / "rarefaction_plots.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        log.info("Saved rarefaction_plots.png")
        for p in existing:
            p.unlink()
            log.info("Removed %s", p.name)
 
    log.info("Done. Output in: %s", out_dir)


if __name__ == "__main__":
    main()
