#!/usr/bin/env python3
"""
select_rarefaction_depth.py

Script to automatically pick a rarefaction depth for microbiome ASV tables.
Has 3 methods: percentile, knee, coverage

Usage:
    python select_rarefaction_depth.py --input asv_table.tsv --output rarefaction_report/ --method knee

Input should be a TSV where rows = samples, columns = taxa.
If it's the other way around it gets transposed automatically.

"""

import argparse
import json
import os
import sys
import warnings
from datetime import datetime
from pathlib import Path
import logging

import numpy as np
import pandas as pd

# setup logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# try to import matplotlib, plots are optional
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    from matplotlib.lines import Line2D
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False
    log.warning("matplotlib not found, skipping plots")


def load_asv_table(path, sep="\t"):
    # load the table
    df = pd.read_csv(path, sep=sep, index_col=0)
    
    # if there are more columns than rows it's probably taxa-as-rows, so transpose
    if df.shape[1] > df.shape[0]:
        log.info("looks like taxa are rows, transposing...")
        df = df.T
    
    df = df.fillna(0).astype(int)
    return df


def goods_coverage(counts):
    """
    Calculates Good's coverage for a sample.
    Formula: 1 - (number of singletons / total reads)
    """
    counts = counts[counts > 0]
    total = counts.sum()
    
    if total == 0:
        return np.nan
    
    singletons = (counts == 1).sum()
    coverage = 1.0 - (singletons / total)
    return coverage


def rarefy_once(counts, depth, rng):
    """subsample reads down to 'depth' without replacement"""
    # expand counts into individual reads, then sample
    reads = np.repeat(np.arange(len(counts)), counts)
    chosen = rng.choice(reads, size=depth, replace=False)
    result = np.bincount(chosen, minlength=len(counts))
    return result


def observed_richness(counts, depth, iterations, rng):
    """
    rarefies the sample multiple times and returns mean richness
    (number of ASVs observed)
    """
    all_richness = []
    for i in range(iterations):
        rarefied = rarefy_once(counts, depth, rng)
        n_asvs = (rarefied > 0).sum()
        all_richness.append(n_asvs)
    
    return float(np.mean(all_richness))


def find_knee(values):
    """
    Find the knee point in a sorted array using the chord/distance method.
    Basically finds the point furthest from the line connecting the first and last points.
    """
    sorted_vals = np.sort(values)
    n = len(sorted_vals)
    
    if n < 3:
        # not enough points
        return float(np.percentile(values, 10))
    
    # normalize to 0-1 range
    xs = np.linspace(0, 1, n)
    y_min = sorted_vals.min()
    y_max = sorted_vals.max()
    
    if y_max == y_min:
        return float(y_min)
    
    ys = (sorted_vals - y_min) / (y_max - y_min)
    
    # direction vector of chord from first to last point
    dx = xs[-1] - xs[0]
    dy = ys[-1] - ys[0]
    
    # perpendicular distance from each point to the chord line
    # formula from analytical geometry
    distances = np.abs(dy * xs - dx * ys + xs[-1] * ys[0] - ys[-1] * xs[0]) / (np.sqrt(dx**2 + dy**2) + 1e-12)
    
    knee_idx = int(np.argmax(distances))
    return float(sorted_vals[knee_idx])


# threshold selection methods

def select_percentile(lib_sizes, percentile):
    threshold = int(np.percentile(lib_sizes, percentile))
    log.info("percentile threshold (%.0f%%): %d", percentile, threshold)
    return threshold


def select_knee(lib_sizes):
    threshold = int(find_knee(lib_sizes))
    log.info("knee threshold: %d", threshold)
    return threshold


def select_coverage(asv_table, lib_sizes, coverage_target, dropout_max, step):
    """
    Find the lowest depth where enough samples still hit the coverage target.
    Steps through possible depths and checks what fraction pass.
    """
    rng = np.random.seed(18)
    max_depth = lib_sizes.max()
    depths = np.arange(step, max_depth + step, step)
    
    passing_fracs = []
    
    for d in depths:
        # only consider samples that have enough reads for this depth
        eligible_mask = lib_sizes >= d
        eligible = asv_table[eligible_mask]
        
        if len(eligible) == 0:
            passing_fracs.append(0.0)
            continue
        
        coverages = []
        for _, row in eligible.iterrows():
            counts = row.values
            rarefied = rarefy_once(counts, int(d), rng)
            cov = goods_coverage(rarefied)
            coverages.append(cov)
        
        # fraction of samples that pass the coverage target
        valid_covs = [c for c in coverages if not np.isnan(c)]
        frac = np.mean([c >= coverage_target for c in valid_covs])
        passing_fracs.append(frac)
    
    passing_fracs = np.array(passing_fracs)
    good_depths = depths[passing_fracs >= (1 - dropout_max)]
    
    if len(good_depths) == 0:
        log.warning("coverage target never met, falling back to 10th percentile")
        return select_percentile(lib_sizes, 10)
    
    threshold = int(good_depths.min())
    log.info("coverage threshold: %d", threshold)
    return threshold


def compute_curves(asv_table, lib_sizes, pass_mask, step, iterations):
    """
    Compute rarefaction curves for all samples.
    Returns a dataframe with columns: sample, depth, richness, passes
    """
    rng = np.random.seed(18)
    records = []
    
    for i, (sample_name, row) in enumerate(asv_table.iterrows()):
        counts = row.values
        N = lib_sizes[i]
        
        # build list of depths to evaluate for this sample
        depths = list(range(step, int(N), step)) + [int(N)]
        
        for d in depths:
            # cap iterations for speed
            richness = observed_richness(counts, d, min(iterations, 20), rng)
            records.append({
                "sample": sample_name,
                "depth": d,
                "richness": richness,
                "passes": bool(pass_mask[i])
            })
    
    return pd.DataFrame(records)


# plotting functions

def plot_curves(curve_df, threshold, out_dir):
    fig, ax = plt.subplots(figsize=(10, 5))
    
    for sample_name, grp in curve_df.groupby("sample"):
        if grp["passes"].iloc[0]:
            colour = "#2980B9"
            alpha = 0.6
        else:
            colour = "#BDC3C7"
            alpha = 0.25
        ax.plot(grp["depth"], grp["richness"], color=colour, alpha=alpha, linewidth=0.7)
    
    ax.axvline(threshold, color="#E74C3C", linestyle="--", linewidth=1.2)
    ax.set_xlabel("Sequencing depth")
    ax.set_ylabel("Observed ASVs")
    ax.set_title("Rarefaction curves per sample")
    
    legend_handles = [
        Line2D([0], [0], color="#2980B9", linewidth=1.5, label="Retained"),
        Line2D([0], [0], color="#BDC3C7", linewidth=1.5, label="Excluded"),
        Line2D([0], [0], color="#E74C3C", linestyle="--", linewidth=1.5, label=f"threshold = {threshold:,}"),
    ]
    ax.legend(handles=legend_handles, frameon=False)
    fig.tight_layout()
    fig.savefig(out_dir / "rarefaction_curves.pdf", dpi=150)
    plt.close(fig)


def plot_libsize(lib_sizes, pass_mask, threshold, out_dir):
    fig, ax = plt.subplots(figsize=(7, 4))
    
    ax.hist(lib_sizes[pass_mask],  bins=40, color="#2980B9", alpha=0.7, label="Retained")
    ax.hist(lib_sizes[~pass_mask], bins=40, color="#BDC3C7", alpha=0.7, label="Excluded")
    ax.axvline(threshold, color="#E74C3C", linestyle="--", linewidth=1.2)
    
    ax.set_xlabel("Read count")
    ax.set_ylabel("Number of samples")
    ax.set_title("Library size distribution")
    ax.legend(frameon=False)
    
    # format x axis with commas
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    
    fig.tight_layout()
    fig.savefig(out_dir / "library_size_distribution.pdf", dpi=150)
    plt.close(fig)


def plot_coverage(qc_df, threshold, out_dir):
    # only look at samples that passed
    cov_vals = qc_df.loc[qc_df["passes_threshold"], "coverage"].dropna()
    
    if len(cov_vals) == 0:
        return
    
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(cov_vals, bins=30, color="#27AE60", alpha=0.8)
    ax.axvline(0.99, color="#E74C3C", linestyle="--", linewidth=1.2, label="99% target")
    
    ax.set_xlabel("Good's coverage")
    ax.set_ylabel("Number of samples")
    ax.set_title(f"Good's coverage at depth = {threshold:,}")
    ax.legend(frameon=False)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: f"{x:.0%}"))
    
    fig.tight_layout()
    fig.savefig(out_dir / "coverage_at_threshold.pdf", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Pick a rarefaction depth automatically from an ASV table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True,  help="ASV table (.tsv or .csv)")
    parser.add_argument("--output", default="rarefaction_report", help="Output folder")
    parser.add_argument("--method", default="knee", choices=["percentile", "knee", "coverage"])
    parser.add_argument("--percentile", type=float, default=10,   help="Which percentile to use (for percentile method)")
    parser.add_argument("--coverage-target", type=float, default=0.99, dest="coverage_target", help="Minimum Good's coverage (for coverage method)")
    parser.add_argument("--dropout-max", type=float, default=0.10, dest="dropout_max",     help="Max fraction of samples to drop")
    parser.add_argument("--step", type=int, default=500, help="Step size when building rarefaction curves")
    parser.add_argument("--iterations", type=int, default=100, help="How many times to rarefy per depth (for curve smoothing)")
    parser.add_argument("--sep", default="\t", help="Separator character (tab or comma)")
    
    args = parser.parse_args()

    # make output dir
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    # step 1: load table
    log.info("Step 1/5: loading table: %s", args.input)
    asv = load_asv_table(args.input, sep=args.sep)
    
    lib_sizes = asv.sum(axis=1).values
    n_samples = len(asv)
    
    log.info("  samples : %d", n_samples)
    log.info("  min lib : %d", lib_sizes.min())
    log.info("  max lib : %d", lib_sizes.max())
    log.info("  median  : %d", int(np.median(lib_sizes)))

    # step 2: pick threshold
    log.info("Step 2/5: selecting threshold (method=%s)", args.method)
    
    if args.method == "percentile":
        threshold = select_percentile(lib_sizes, args.percentile)
    elif args.method == "knee":
        threshold = select_knee(lib_sizes)
    else:
        threshold = select_coverage(asv, lib_sizes, args.coverage_target, args.dropout_max, args.step)

    pass_mask = lib_sizes >= threshold
    n_pass = int(pass_mask.sum())
    n_fail = n_samples - n_pass
    dropout_frac = n_fail / n_samples

    log.info("  kept %d / %d samples (dropped %.1f%%)", n_pass, n_samples, dropout_frac * 100)

    if dropout_frac > args.dropout_max:
        warnings.warn(
            f"dropout rate ({dropout_frac:.1%}) is higher than --dropout-max ({args.dropout_max:.1%}). "
            "You might want to lower the threshold or remove low-depth samples manually."
        )

    # step 3: rarefaction curves
    log.info("Step 3/5: computing rarefaction curves...")
    curve_df = compute_curves(asv, lib_sizes, pass_mask, args.step, args.iterations)

    # step 4: plots
    if HAS_PLOT:
        log.info("Step 4/5: generating plots...")
        plot_curves(curve_df, threshold, out_dir)
        plot_libsize(lib_sizes, pass_mask, threshold, out_dir)
    else:
        log.info("Step 4/5: skipping plots (matplotlib not available)")

    # compute per-sample coverage at the chosen threshold
    rng = np.random.seed(18)
    coverages = []
    
    for i, (_, row) in enumerate(asv.iterrows()):
        counts = row.values
        if lib_sizes[i] >= threshold:
            rarefied = rarefy_once(counts, threshold, rng)
            cov = goods_coverage(rarefied)
            coverages.append(cov)
        else:
            coverages.append(float("nan"))  # can't rarefy below threshold

    qc_df = pd.DataFrame({
        "sample":           asv.index,
        "lib_size":         lib_sizes,
        "passes_threshold": pass_mask,
        "coverage":         coverages,
    })

    if HAS_PLOT:
        plot_coverage(qc_df, threshold, out_dir)

    # step 5: save outputs
    log.info("Step 5/5: writing outputs to %s", out_dir)
    
    (out_dir / "rarefaction_threshold.txt").write_text(str(threshold) + "\n")
    qc_df.to_csv(out_dir / "sample_qc.tsv", sep="\t", index=False)
    curve_df.to_csv(out_dir / "rarefaction_curves.tsv", sep="\t", index=False)

    # summary json for downstream tools
    summary = {
        "date":            datetime.now().isoformat(),
        "input":           args.input,
        "method":          args.method,
        "n_samples":       n_samples,
        "lib_size_min":    int(lib_sizes.min()),
        "lib_size_max":    int(lib_sizes.max()),
        "lib_size_median": int(np.median(lib_sizes)),
        "threshold":       threshold,
        "n_retained":      n_pass,
        "n_dropped":       n_fail,
        "dropout_frac":    round(dropout_frac, 4),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # print a text report
    report_lines = [
        "=== Rarefaction threshold selection report ===",
        f"Date        : {summary['date']}",
        f"Input file  : {args.input}",
        f"Method      : {args.method}",
        "",
        "--- Library size summary ---",
        f"  Samples   : {n_samples}",
        f"  Min       : {lib_sizes.min():,}",
        f"  Max       : {lib_sizes.max():,}",
        f"  Mean      : {int(lib_sizes.mean()):,}",
        f"  Median    : {int(np.median(lib_sizes)):,}",
        f"  10th pct  : {int(np.percentile(lib_sizes, 10)):,}",
        "",
        "--- Selected threshold ---",
        f"  Threshold : {threshold:,} reads",
        f"  Retained  : {n_pass} / {n_samples} samples",
        f"  Dropped   : {n_fail} samples ({dropout_frac:.1%})",
        "",
        "--- Output files ---",
        "  rarefaction_threshold.txt",
        "  sample_qc.tsv",
        "  summary.json",
        "  rarefaction_curves.tsv",
        "  rarefaction_curves.pdf          (if matplotlib available)",
        "  library_size_distribution.pdf",
        "  coverage_at_threshold.pdf",
    ]
    report_text = "\n".join(report_lines)
    (out_dir / "report.txt").write_text(report_text + "\n")
    print(report_text)

    log.info("done!")


if __name__ == "__main__":
    main()
