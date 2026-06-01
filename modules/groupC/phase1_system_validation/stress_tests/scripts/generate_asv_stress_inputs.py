#!/usr/bin/env python3
"""
generate_asv_stress_inputs.py
─────────────────────────────
Generates stress-test datasets from the PRJEB10949 mock community dataset
(BEI mock communities, Lluch et al. 2015), whose composition is fully known.

Output structure:
    stress_inputs/
    ├── ST00_baseline/                  unmodified asv_table.tsv
    ├── ST01_zero_count_sample/         one even sample forced to 0 reads
    ├── ST02_low_depth/                 biological samples downsampled to 100 reads
    ├── ST04_single_taxon/              even samples collapsed to one dominant ASV
    ├── ST06a_contamination_blanks/     Pseudomonas spiked at 50% in blanks (for decontam)
    └── ST06b_contamination_biological/ Pseudomonas spiked at 30% in even/staggered (for F1)
    technical_checks/
    ├── ST03_invalid_count_table/       non-numeric values in count table (no F1/recall)
    └── ST05_metadata_mismatch/         wrong sample IDs in metadata (no F1/recall)
    ASV_taxonomy.tsv                    shared taxonomy file (unchanged across scenarios)
    metadata_PRJEB10949.tsv             sample metadata with clear columns

Usage:
    python generate_asv_stress_inputs.py \
        --asv_table  results/pipeline_outputs/asv_table.tsv \
        --taxonomy   results/pipeline_outputs/ASV_taxonomy.tsv \
        --outdir     stress_tests/
"""

import argparse
import shutil
from pathlib import Path
import numpy as np
import pandas as pd

SEED = 42
rng  = np.random.default_rng(SEED)

# Sample IDs from PRJEB10949 mock community
EVEN_SAMPLES      = ["ERR1049996", "ERR1049997", "ERR1049998"]
STAGGERED_SAMPLES = ["ERR1049999", "ERR1050000", "ERR1050001"]
BLANK_SAMPLES     = ["ERR1049992", "ERR1049993", "ERR1049994",
                     "ERR1049995", "ERR1049938", "ERR1049939", "ERR1049940"]


def load_asv(path):
    asv = pd.read_csv(path, sep="\t", index_col=0)
    asv.index.name = "#OTU ID"
    return asv


def save_asv(df, path):
    df.index.name = "#OTU ID"
    df.to_csv(path, sep="\t")


def make_dir(outdir, sid, name):
    d = outdir / f"{sid}_{name}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def build_metadata():
    """Build improved metadata with clear columns for all 13 PRJEB10949 samples."""
    rows = []
    for i, s in enumerate(EVEN_SAMPLES):
        rows.append({"sample_id": s, "sample_type": "biological",
                     "community_type": "even", "replicate": i + 1,
                     "is_negative_control": "no"})
    for i, s in enumerate(STAGGERED_SAMPLES):
        rows.append({"sample_id": s, "sample_type": "biological",
                     "community_type": "staggered", "replicate": i + 1,
                     "is_negative_control": "no"})
    for i, s in enumerate(BLANK_SAMPLES):
        rows.append({"sample_id": s, "sample_type": "blank",
                     "community_type": "negative_control", "replicate": i + 1,
                     "is_negative_control": "yes"})
    return pd.DataFrame(rows)


# ── Main stress scenarios (produce F1 / recall / RMSE) ────────────────────────

def ST00_baseline(asv, outdir):
    """Positive control — unmodified data. Reference for all other scenarios."""
    d = make_dir(outdir, "ST00", "baseline")
    save_asv(asv, d / "asv_table.tsv")
    print(f"  ST00 ✓  {asv.shape[1]} samples x {asv.shape[0]} ASVs — unmodified")


def ST01_zero_count_sample(asv, outdir):
    """One even mock replicate forced to 0 reads."""
    d      = make_dir(outdir, "ST01", "zero_count_sample")
    mod    = asv.copy()
    target = EVEN_SAMPLES[0]
    orig   = int(mod[target].sum())
    mod[target] = 0
    save_asv(mod, d / "asv_table.tsv")
    print(f"  ST01 ✓  {target}: {orig} → 0 reads")


def ST02_low_depth(asv, outdir, target_reads=100):
    """All biological mock samples downsampled to 100 reads (below threshold of 316)."""
    d   = make_dir(outdir, "ST02", "low_depth")
    mod = asv.copy()
    for s in EVEN_SAMPLES + STAGGERED_SAMPLES:
        total = int(mod[s].sum())
        if total > target_reads:
            probs  = mod[s].values / total
            mod[s] = rng.multinomial(target_reads, probs)
    save_asv(mod, d / "asv_table.tsv")
    print(f"  ST02 ✓  6 biological samples → {target_reads} reads "
          f"(rarefaction threshold: 316)")


def ST04_single_taxon(asv, outdir):
    """All reads in even samples collapsed into the single most abundant ASV."""
    d   = make_dir(outdir, "ST04", "single_taxon")
    mod = asv.copy()
    for s in EVEN_SAMPLES:
        total  = int(mod[s].sum())
        winner = mod[s].idxmax()
        mod[s] = 0
        mod.loc[winner, s] = total
    save_asv(mod, d / "asv_table.tsv")
    print(f"  ST04 ✓  even samples: all reads → {mod[EVEN_SAMPLES[0]].idxmax()}")


def get_contaminant_asv(taxonomy_path, asv):
    """Find the Pseudomonas ASV in the taxonomy table."""
    tax    = pd.read_csv(taxonomy_path, sep="\t")
    pseudo = tax[tax["Genus"].str.contains("Pseudomonas", case=False, na=False)]
    return pseudo.iloc[0]["ASV_ID"] if len(pseudo) > 0 else asv.index[0]


def ST06a_contamination_blanks(asv, taxonomy_path, outdir, contam_pct=0.5):
    """
    Pseudomonas spiked at 50% ONLY in the H2O negative controls (blanks).
    Purpose: validate Contamination filtering module
             (contamination_filtering.R / Kraken2).
    Note: does NOT affect precision/recall of even/staggered because
          validation_metrics.py does not evaluate blanks.
    """
    d          = make_dir(outdir, "ST06a", "contamination_blanks")
    mod        = asv.copy()
    contam_asv = get_contaminant_asv(taxonomy_path, asv)

    for blank in BLANK_SAMPLES:
        total = int(mod[blank].sum())
        spike = int(total * contam_pct) if total > 0 else 5000
        mod[blank] = (mod[blank] // 2).astype(int)
        mod.loc[contam_asv, blank] += spike

    save_asv(mod, d / "asv_table.tsv")
    pd.DataFrame([{"ASV_ID": contam_asv, "spike_pct": contam_pct,
                   "spike_target": "blanks_only",
                   "targets": ", ".join(BLANK_SAMPLES)}
                  ]).to_csv(d / "spiked_contaminant_ASVs.tsv", sep="\t", index=False)

    rel = mod.loc[contam_asv, BLANK_SAMPLES[0]] / mod[BLANK_SAMPLES[0]].sum()
    print(f"  ST06a ✓  {contam_asv} (Pseudomonas) at {rel:.0%} in blanks "
          f"→ for Pau Villen (decontam)")


def ST06b_contamination_biological(asv, taxonomy_path, outdir, contam_pct=0.3):
    """
    Pseudomonas spiked at 30% in biological mock samples (even + staggered).
    Purpose: measure how contamination affects precision/recall/F1 using
             validation_metrics.py.
    Pseudomonas is not in the even/staggered ground truth, so it appears as
    a false positive → lowers precision.
    """
    d          = make_dir(outdir, "ST06b", "contamination_biological")
    mod        = asv.copy()
    contam_asv = get_contaminant_asv(taxonomy_path, asv)

    for sample in EVEN_SAMPLES + STAGGERED_SAMPLES:
        total = int(mod[sample].sum())
        spike = int(total * contam_pct)
        mod[sample] = (mod[sample] // (1 + contam_pct)).astype(int)
        mod.loc[contam_asv, sample] += spike

    save_asv(mod, d / "asv_table.tsv")
    pd.DataFrame([{"ASV_ID": contam_asv, "spike_pct": contam_pct,
                   "spike_target": "biological_samples",
                   "targets": ", ".join(EVEN_SAMPLES + STAGGERED_SAMPLES)}
                  ]).to_csv(d / "spiked_contaminant_ASVs.tsv", sep="\t", index=False)

    rel = mod.loc[contam_asv, EVEN_SAMPLES[0]] / mod[EVEN_SAMPLES[0]].sum()
    print(f"  ST06b ✓  {contam_asv} (Pseudomonas) at {rel:.0%} in even/staggered "
          f"→ for validation_metrics.py (lowers precision)")


# ── Technical checks (do not produce F1 / recall) ────────────────────

def ST03_invalid_count_table(asv, outdir):
    """10 cells set to NOT_A_NUMBER — tests whether the pipeline detects corrupt input."""
    d   = make_dir(outdir, "ST03", "invalid_count_table")
    mod = asv.copy().astype(object)
    for row in rng.choice(mod.index, size=10, replace=False):
        mod.loc[row, STAGGERED_SAMPLES[0]] = "NOT_A_NUMBER"
    save_asv(mod, d / "asv_table.tsv")
    print(f"  ST03 ✓  10 cells → 'NOT_A_NUMBER' (technical check)")


def ST05_metadata_mismatch(asv, outdir):
    """Sample IDs in metadata have a _WRONG suffix — tests mismatch detection."""
    d    = make_dir(outdir, "ST05", "metadata_mismatch")
    save_asv(asv, d / "asv_table.tsv")
    meta = build_metadata().copy()
    meta["sample_id"] = meta["sample_id"] + "_WRONG"
    meta.to_csv(d / "metadata_wrong.tsv", sep="\t", index=False)
    print(f"  ST05 ✓  metadata IDs with '_WRONG' suffix (0/13 match) (technical check)")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--asv_table", required=True,
                    help="Path to results/pipeline_outputs/asv_table.tsv")
    ap.add_argument("--taxonomy",  required=True,
                    help="Path to results/pipeline_outputs/ASV_taxonomy.tsv")
    ap.add_argument("--outdir",    required=True,
                    help="Output directory for stress inputs")
    args = ap.parse_args()

    outdir   = Path(args.outdir)
    main_dir = outdir / "stress_inputs"
    tech_dir = outdir / "optional_technical_checks"
    main_dir.mkdir(parents=True, exist_ok=True)
    tech_dir.mkdir(parents=True, exist_ok=True)

    asv = load_asv(args.asv_table)

    print("=" * 60)
    print("  Stress Test Dataset Generator — Group C / Issue #9")
    print("  Dataset: PRJEB10949 (BEI mock communities, Lluch 2015)")
    print("=" * 60)
    print(f"\n  Input: {asv.shape[1]} samples x {asv.shape[0]} ASVs\n")

    # Shared files written once to the root output directory
    shutil.copy(args.taxonomy, outdir / "ASV_taxonomy.tsv")
    build_metadata().to_csv(outdir / "metadata_PRJEB10949.tsv", sep="\t", index=False)
    print(f"  taxonomy ✓  ASV_taxonomy.tsv")
    print(f"  metadata ✓  metadata_PRJEB10949.tsv\n")

    print("  -- Main scenarios (produce F1 / recall / RMSE) --")
    ST00_baseline(asv, main_dir)
    ST01_zero_count_sample(asv, main_dir)
    ST02_low_depth(asv, main_dir)
    ST04_single_taxon(asv, main_dir)
    ST06a_contamination_blanks(asv, args.taxonomy, main_dir)
    ST06b_contamination_biological(asv, args.taxonomy, main_dir)

    print("\n  -- Optional technical checks --")
    ST03_invalid_count_table(asv, tech_dir)
    ST05_metadata_mismatch(asv, tech_dir)

    print(f"\n  ✓ stress_inputs/             → {main_dir}")
    print(f"  ✓ optional_technical_checks/ → {tech_dir}")
    print(f"  ✓ ASV_taxonomy.tsv")
    print(f"  ✓ metadata_PRJEB10949.tsv")
    print(f"\n  Next step: run validation_metrics.py on each scenario")


if __name__ == "__main__":
    main()
