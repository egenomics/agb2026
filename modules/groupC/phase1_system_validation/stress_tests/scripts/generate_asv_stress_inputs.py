#!/usr/bin/env python3
"""
Generate small ASV-table stress-test inputs from Group B output files.

Usage:
    python scripts/generate_asv_stress_inputs.py \
        --asv_table asv_table.tsv \
        --taxonomy ASV_taxonomy.tsv \
        --outdir stress_inputs_regenerated

This script creates small modified ASV tables for Group C robustness testing:
- valid subset
- zero-count sample
- very low-depth samples
- invalid non-numeric count
- single-taxon dominance
- metadata mismatch
- contamination spike
"""

from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import shutil

def write_metadata(outpath, sample_ids, test_id, sample_type="sample", extra_notes=""):
    if isinstance(sample_type, str):
        sample_type = [sample_type] * len(sample_ids)
    if isinstance(extra_notes, str):
        extra_notes = [extra_notes] * len(sample_ids)

    md = pd.DataFrame({
        "sample-id": sample_ids,
        "stress_test_id": test_id,
        "sample_type": sample_type,
        "is_negative_control": ["yes" if x in ["negative_control", "blank"] else "no" for x in sample_type],
        "notes": extra_notes
    })
    md.to_csv(outpath, sep="\t", index=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asv_table", required=True)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)

    rng = np.random.default_rng(args.seed)
    asv = pd.read_csv(args.asv_table, sep="\t")
    tax = pd.read_csv(args.taxonomy, sep="\t")
    id_col = asv.columns[0]
    sample_cols = list(asv.columns[1:])
    counts = asv[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    sample_totals = counts.sum(axis=0).sort_values()

    top_asv_ids = counts.sum(axis=1).sort_values(ascending=False).head(50).index
    valid_samples = sample_totals[sample_totals > 10000].index[:10].tolist()
    subset = pd.concat([asv.loc[top_asv_ids, [id_col]], counts.loc[top_asv_ids, valid_samples]], axis=1)
    subset_tax = tax[tax["ASV_ID"].isin(subset[id_col])].copy()

    # ST00
    d = outdir / "ST00_valid_subset"
    d.mkdir()
    subset.to_csv(d / "asv_table_subset.tsv", sep="\t", index=False)
    subset_tax.to_csv(d / "ASV_taxonomy_subset.tsv", sep="\t", index=False)
    write_metadata(d / "metadata_valid_subset.tsv", valid_samples, "ST00", "sample", "Valid subset control")

    # ST01
    d = outdir / "ST01_zero_count_sample"
    d.mkdir()
    zero_table = subset.copy()
    zero_table["ST01_zero_reads"] = 0
    zero_cols = valid_samples[:5] + ["ST01_zero_reads"]
    zero_table = zero_table[[id_col] + zero_cols]
    zero_table.to_csv(d / "asv_table_zero_count_sample.tsv", sep="\t", index=False)
    write_metadata(d / "metadata_zero_count_sample.tsv", zero_cols, "ST01", "sample", "Includes one zero-depth sample")

    # ST02
    d = outdir / "ST02_low_depth"
    d.mkdir()
    low_cols = valid_samples[:5]
    low_table = pd.DataFrame({id_col: subset[id_col]})
    for col in low_cols:
        v = subset[col].to_numpy(dtype=int)
        probs = v / v.sum()
        low_table[f"{col}_downsampled_100"] = rng.multinomial(100, probs)
    low_table.to_csv(d / "asv_table_low_depth_100_reads.tsv", sep="\t", index=False)
    write_metadata(d / "metadata_low_depth_100_reads.tsv", list(low_table.columns[1:]), "ST02", "sample", "Artificially downsampled to 100 reads")

    # ST03
    d = outdir / "ST03_invalid_count_table"
    d.mkdir()
    bad_table = subset[[id_col] + valid_samples[:3]].copy()
    bad_table.iloc[0, 1] = "NOT_A_NUMBER"
    bad_table.to_csv(d / "asv_table_invalid_non_numeric_count.tsv", sep="\t", index=False)
    write_metadata(d / "metadata_invalid_count_table.tsv", valid_samples[:3], "ST03", "sample", "One non-numeric count value")

    # ST04
    d = outdir / "ST04_single_taxon"
    d.mkdir()
    tax_counts = tax.merge(pd.DataFrame({"ASV_ID": asv[id_col], "total": counts.sum(axis=1).values}), on="ASV_ID")
    bact = tax_counts[tax_counts["Genus"].fillna("").str.contains("Bacteroides", case=False, na=False)].sort_values("total", ascending=False)
    single_asv = bact.iloc[0]["ASV_ID"] if len(bact) else tax_counts.sort_values("total", ascending=False).iloc[0]["ASV_ID"]
    single_table = subset[[id_col]].copy()
    single_cols = [f"ST04_single_taxon_{i}" for i in range(1, 6)]
    for c in single_cols:
        single_table[c] = 0
    single_table.loc[single_table[id_col] == single_asv, single_cols] = 10000
    single_table.to_csv(d / "asv_table_single_taxon_dominance.tsv", sep="\t", index=False)
    subset_tax.to_csv(d / "ASV_taxonomy_single_taxon_subset.tsv", sep="\t", index=False)
    write_metadata(d / "metadata_single_taxon_dominance.tsv", single_cols, "ST04", "sample", f"All reads assigned to {single_asv}")

    # ST05
    d = outdir / "ST05_metadata_mismatch"
    d.mkdir()
    subset[[id_col] + valid_samples[:5]].to_csv(d / "asv_table_for_metadata_mismatch.tsv", sep="\t", index=False)
    mismatch_ids = [f"{s}_WRONG_ID" for s in valid_samples[:5]]
    write_metadata(d / "metadata_mismatch_wrong_sample_ids.tsv", mismatch_ids, "ST05", "sample", "Metadata IDs do not match ASV table")

    # ST06
    d = outdir / "ST06_contamination_spike"
    d.mkdir()
    contam_table = subset[[id_col] + valid_samples[:5]].copy()
    candidate_genera = ["Pseudomonas", "Escherichia-Shigella", "Klebsiella", "Streptococcus"]
    contam_asvs = []
    for gen in candidate_genera:
        sub = tax_counts[tax_counts["Genus"].fillna("").str.fullmatch(gen, case=False, na=False)].sort_values("total", ascending=False)
        if len(sub):
            contam_asvs.append(sub.iloc[0]["ASV_ID"])
    contam_asvs = list(dict.fromkeys(contam_asvs))[:3]
    if not contam_asvs:
        contam_asvs = tax_counts.sort_values("total", ascending=False).head(3)["ASV_ID"].tolist()

    for nc in ["NC01", "NC02", "NC03"]:
        contam_table[nc] = 0
    for asv_id in contam_asvs:
        contam_table.loc[contam_table[id_col] == asv_id, ["NC01", "NC02", "NC03"]] = [5000, 7000, 6000]
    contam_table.to_csv(d / "asv_table_contamination_spike.tsv", sep="\t", index=False)
    subset_tax.to_csv(d / "ASV_taxonomy_contamination_spike_subset.tsv", sep="\t", index=False)
    sample_ids = valid_samples[:5] + ["NC01", "NC02", "NC03"]
    sample_types = ["sample"] * 5 + ["negative_control"] * 3
    notes = ["Biological sample"] * 5 + [f"Negative control spiked with {','.join(contam_asvs)}"] * 3
    write_metadata(d / "metadata_contamination_spike.tsv", sample_ids, "ST06", sample_types, notes)
    tax[tax["ASV_ID"].isin(contam_asvs)].to_csv(d / "spiked_contaminant_ASVs.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
