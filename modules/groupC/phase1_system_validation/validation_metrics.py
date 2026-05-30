#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from scipy.spatial.distance import braycurtis


# ===============
# PARSE ARGUMENTS
# ===============
if len(sys.argv) != 4:
    print("Usage: python validation_metrics.py <ground_truth> <asv_table> <taxonomy>")
    sys.exit(1)

ground_truth_file = sys.argv[1]
asv_table_file    = sys.argv[2]
taxonomy_file     = sys.argv[3]


# ==================
# SAMPLE DEFINITIONS
# ==================
even_replicates      = ["ERR1049996", "ERR1049997", "ERR1049998"]
staggered_replicates = ["ERR1049999", "ERR1050000", "ERR1050001"]
all_replicates       = even_replicates + staggered_replicates

# Species excluded from FN calculation
excluded_from_fn = ["Deinococcus radiodurans"]
excluded_genera  = [s.split()[0] for s in excluded_from_fn]


# =========
# LOAD DATA
# =========
ground_truth = pd.read_csv(ground_truth_file, sep="\t")
asv_table    = pd.read_csv(asv_table_file, sep="\t")
asv_table    = asv_table.rename(columns={"#OTU ID": "ASV_ID"})

# Handle comment="#" not catching "#OTU ID" header
if asv_table.columns[0] != "ASV_ID":
    asv_table = asv_table.rename(columns={asv_table.columns[0]: "ASV_ID"})

taxonomy = pd.read_csv(taxonomy_file, sep="\t")


# ===========================
# EXTRACT GROUND TRUTH GENERA
# ===========================
# Output from group B is not at species level --> genus
ground_truth["genus"] = ground_truth["species"].str.split().str[0]

gt_species = set(ground_truth["species"]) - set(excluded_from_fn)
gt_genera  = set(ground_truth["genus"])   - set(excluded_genera)


# =============================
# MERGE TAXONOMY INTO ASV TABLE
# =============================
# Genus values like "Escherichia-Shigella" — take only the first part before hyphen
taxonomy["Genus_clean"] = taxonomy["Genus"].str.split("-").str[0]

merged = asv_table.merge(
    taxonomy[["ASV_ID", "Genus", "Genus_clean"]],
    on="ASV_ID",
    how="left"
)


# ================
# HELPER FUNCTIONS
# ================
def compute_metrics(tp, fp, fn):
    precision = tp / (tp + fp)           if (tp + fp) > 0           else 0.0
    recall    = tp / (tp + fn)           if (tp + fn) > 0           else 0.0
    f1        = (2 * precision * recall) / (precision + recall)     if (precision + recall) > 0 else 0.0
    accuracy  = tp / (tp + fp + fn)      if (tp + fp + fn) > 0      else 0.0
    return round(precision, 4), round(recall, 4), round(f1, 4), round(accuracy, 4)

def normalize(series):
    total = series.sum()
    return series / total if total > 0 else series

def compute_rmse(expected, observed):
    return round(np.sqrt(np.mean((expected - observed) ** 2)), 4)

def compute_bray_curtis(expected, observed):
    return round(braycurtis(expected, observed), 4)


# ===============
# COMPUTE METRICS
# ===============
detection_rows = []
abundance_rows = []

for replicate in all_replicates:

    if replicate not in merged.columns:
        print(f"Warning: {replicate} not found in ASV table, skipping.")
        continue

    community_type = "even" if replicate in even_replicates else "staggered"
    pct_col        = "even_expected_pct" if community_type == "even" else "staggered_expected_pct"

    # --- GENUS LEVEL ---
    detected_genera = set(
        merged.loc[merged[replicate] > 0, "Genus_clean"].dropna()
    )

    tp_g = len(detected_genera & gt_genera)
    fp_g = len(detected_genera - gt_genera)
    fn_g = len(gt_genera - detected_genera)
    precision_g, recall_g, f1_g, accuracy_g = compute_metrics(tp_g, fp_g, fn_g)

    detection_rows.append({
        "replicate"      : replicate,
        "community_type" : community_type,
        "level"          : "genus",
        "TP"             : tp_g,
        "FP"             : fp_g,
        "FN"             : fn_g,
        "accuracy"       : accuracy_g,
        "precision"      : precision_g,
        "recall"         : recall_g,
        "f1"             : f1_g
    })

    # --- ABUNDANCE METRICS ---
    # Aggregate counts by genus, normalize to relative abundance
    observed_counts = merged[[replicate, "Genus_clean"]].copy()
    observed_counts = observed_counts.groupby("Genus_clean")[replicate].sum()
    observed_rel    = normalize(observed_counts)

    # Align with ground truth at genus level
    gt_abundance = ground_truth.groupby("genus")[pct_col].sum() / 100
    all_genera   = gt_abundance.index.union(observed_rel.index)
    gt_aligned   = gt_abundance.reindex(all_genera, fill_value=0)
    obs_aligned  = observed_rel.reindex(all_genera, fill_value=0)

    rmse = compute_rmse(gt_aligned.values, obs_aligned.values)
    bc   = compute_bray_curtis(gt_aligned.values, obs_aligned.values)

    abundance_rows.append({
        "replicate"      : replicate,
        "community_type" : community_type,
        "RMSE"           : rmse,
        "bray_curtis"    : bc
    })


# ================
# COMPUTE AVERAGES
# ================
detection_df = pd.DataFrame(detection_rows)
abundance_df = pd.DataFrame(abundance_rows)

det_avg = detection_df.groupby(["community_type", "level"])[
    ["accuracy", "precision", "recall", "f1"]
].mean().round(4).reset_index()
det_avg["replicate"] = "average"

abd_avg = abundance_df.groupby("community_type")[
    ["RMSE", "bray_curtis"]
].mean().round(4).reset_index()
abd_avg["replicate"] = "average"

detection_df = pd.concat([detection_df, det_avg], ignore_index=True)
abundance_df = pd.concat([abundance_df, abd_avg], ignore_index=True)


# ======
# EXPORT
# ======
detection_df.to_csv("detection_metrics.tsv", sep="\t", index=False)
abundance_df.to_csv("abundance_metrics.tsv",  sep="\t", index=False)
