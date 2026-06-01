#!/usr/bin/env python3
"""
summarize_stress_tests.py

Automated summary for Group C PRJEB10949 stress testing.

Run from:
    modules/groupC/phase1_system_validation/stress_tests

Typical workflow:
    bash scripts/run_validation_metrics_all.sh
    python scripts/summarize_stress_tests.py

Inputs expected:
    stress_results/<scenario>/detection_metrics.tsv
    stress_results/<scenario>/abundance_metrics.tsv
    stress_inputs/ST06a_contamination_blanks/asv_table.tsv
    stress_inputs/ST06a_contamination_blanks/spiked_contaminant_ASVs.tsv
    technical_checks/ST03_invalid_count_table/asv_table.tsv
    technical_checks/ST05_metadata_mismatch/asv_table.tsv
    technical_checks/ST05_metadata_mismatch/metadata_wrong.tsv
    metadata_PRJEB10949.tsv

Outputs:
    stress_test_results.tsv
    technical_check_results.tsv
    stress_test_report.md
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd


MAIN_SCENARIOS = [
    "ST00_baseline",
    "ST01_zero_count_sample",
    "ST02_low_depth",
    "ST04_single_taxon",
    "ST06b_contamination_biological",
]

SCENARIO_DESCRIPTIONS = {
    "ST00_baseline": "Unmodified PRJEB10949 baseline used as reference.",
    "ST01_zero_count_sample": "One even mock replicate was forced to zero reads.",
    "ST02_low_depth": "Biological mock samples were downsampled to low sequencing depth.",
    "ST04_single_taxon": "Even mock samples were collapsed into a single dominant ASV.",
    "ST06b_contamination_biological": "Pseudomonas contamination was spiked into biological mock samples.",
    "ST06a_contamination_blanks": "Pseudomonas contamination was spiked into H2O blank controls.",
    "ST03_invalid_count_table": "Count table contains non-numeric values.",
    "ST05_metadata_mismatch": "Metadata sample IDs do not match ASV table sample IDs.",
}


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_csv(path, sep="\t")


def norm_col(name: str) -> str:
    return name.strip().lower().replace("-", "_").replace(" ", "_").replace(".", "_")


def find_col(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    cols = {norm_col(c): c for c in df.columns}
    for c in candidates:
        nc = norm_col(c)
        if nc in cols:
            return cols[nc]
    return None


def first_col(df: pd.DataFrame) -> str:
    return df.columns[0]


def identify_community_column(df: pd.DataFrame) -> Optional[str]:
    return find_col(df, ["community_type", "community", "group", "mock_type", "sample_type"])


def identify_sample_column(df: pd.DataFrame) -> str:
    return find_col(df, ["sample", "sample_id", "sample-id", "replicate", "run", "run_accession"]) or first_col(df)


def identify_metric_column(df: pd.DataFrame, names: Iterable[str]) -> Optional[str]:
    col = find_col(df, names)
    if col:
        return col
    norm_to_orig = {norm_col(c): c for c in df.columns}
    for target in names:
        nt = norm_col(target)
        for nc, orig in norm_to_orig.items():
            if nt in nc:
                return orig
    return None


def filter_average_rows(df: pd.DataFrame) -> pd.DataFrame:
    sample_col = identify_sample_column(df)
    sample_values = df[sample_col].astype(str).str.lower()
    avg = df[sample_values.str.contains("average", na=False)].copy()
    if not avg.empty:
        return avg

    for c in df.columns:
        vals = df[c].astype(str).str.lower()
        avg = df[vals.str.contains("average", na=False)].copy()
        if not avg.empty:
            return avg

    return df.copy()


def infer_community_from_row(row: pd.Series, community_col: Optional[str], sample_col: str) -> str:
    if community_col and community_col in row.index:
        val = str(row[community_col]).strip().lower()
        if val and val not in {"nan", "none"}:
            if "stagger" in val:
                return "staggered"
            if "even" in val:
                return "even"
            if "blank" in val or "h2o" in val or "negative" in val:
                return "blank"
            return val

    val = str(row[sample_col]).strip().lower()
    if "stagger" in val:
        return "staggered"
    if "even" in val:
        return "even"
    if "blank" in val or "h2o" in val or "negative" in val:
        return "blank"

    return val


def summarise_detection(detection_path: Path, scenario: str) -> pd.DataFrame:
    df = read_tsv(detection_path)

    sample_col = identify_sample_column(df)
    community_col = identify_community_column(df)

    precision_col = identify_metric_column(df, ["precision"])
    recall_col = identify_metric_column(df, ["recall"])
    f1_col = identify_metric_column(df, ["f1", "f1_score", "f1-score"])
    accuracy_col = identify_metric_column(df, ["accuracy", "acc"])

    required = {
        "precision": precision_col,
        "recall": recall_col,
        "f1": f1_col,
        "accuracy": accuracy_col,
    }
    missing = [k for k, v in required.items() if v is None]
    if missing:
        raise ValueError(
            f"{detection_path} is missing expected metric columns: {missing}. "
            f"Found columns: {list(df.columns)}"
        )

    for col in [precision_col, recall_col, f1_col, accuracy_col]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    avg = filter_average_rows(df)
    avg["community_type_inferred"] = avg.apply(
        lambda r: infer_community_from_row(r, community_col, sample_col), axis=1
    )

    has_average_rows = avg[sample_col].astype(str).str.lower().str.contains("average", na=False).any()

    if not has_average_rows:
        grouped = (
            avg.groupby("community_type_inferred", dropna=False)[
                [precision_col, recall_col, f1_col, accuracy_col]
            ]
            .mean(numeric_only=True)
            .reset_index()
        )
    else:
        grouped = avg[["community_type_inferred", precision_col, recall_col, f1_col, accuracy_col]].copy()

    grouped = grouped.rename(
        columns={
            "community_type_inferred": "community_type",
            precision_col: "precision",
            recall_col: "recall",
            f1_col: "f1",
            accuracy_col: "accuracy",
        }
    )
    grouped.insert(0, "test_id", scenario)
    grouped["description"] = SCENARIO_DESCRIPTIONS.get(scenario, "")
    return grouped


def summarise_abundance(abundance_path: Path, scenario: str) -> pd.DataFrame:
    if not abundance_path.exists():
        return pd.DataFrame(columns=["test_id", "community_type", "rmse", "bray_curtis"])

    df = read_tsv(abundance_path)
    sample_col = identify_sample_column(df)
    community_col = identify_community_column(df)

    rmse_col = identify_metric_column(df, ["rmse", "root_mean_square_error"])
    bray_col = identify_metric_column(df, ["bray_curtis", "bray-curtis", "braycurtis"])

    if rmse_col is None and bray_col is None:
        return pd.DataFrame(columns=["test_id", "community_type", "rmse", "bray_curtis"])

    metric_cols = [c for c in [rmse_col, bray_col] if c is not None]
    for col in metric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    avg = filter_average_rows(df)
    avg["community_type_inferred"] = avg.apply(
        lambda r: infer_community_from_row(r, community_col, sample_col), axis=1
    )

    has_average_rows = avg[sample_col].astype(str).str.lower().str.contains("average", na=False).any()

    if not has_average_rows:
        grouped = (
            avg.groupby("community_type_inferred", dropna=False)[metric_cols]
            .mean(numeric_only=True)
            .reset_index()
        )
    else:
        grouped = avg[["community_type_inferred"] + metric_cols].copy()

    grouped = grouped.rename(columns={"community_type_inferred": "community_type"})
    if rmse_col:
        grouped = grouped.rename(columns={rmse_col: "rmse"})
    else:
        grouped["rmse"] = pd.NA
    if bray_col:
        grouped = grouped.rename(columns={bray_col: "bray_curtis"})
    else:
        grouped["bray_curtis"] = pd.NA

    grouped.insert(0, "test_id", scenario)
    return grouped[["test_id", "community_type", "rmse", "bray_curtis"]]


def interpret_stress_row(row: pd.Series) -> str:
    test_id = row.get("test_id", "")
    if test_id == "ST00_baseline":
        return "Baseline reference."

    delta_f1 = row.get("delta_f1_vs_baseline", pd.NA)
    if pd.isna(delta_f1):
        return "No baseline comparison available."

    if delta_f1 <= -0.30:
        return "Strong performance degradation compared with baseline."
    if delta_f1 <= -0.10:
        return "Moderate performance degradation compared with baseline."
    if delta_f1 < -0.02:
        return "Slight performance decrease compared with baseline."
    if delta_f1 <= 0.02:
        return "Performance remained broadly stable compared with baseline."
    return "Performance increased relative to baseline; inspect precision/recall trade-off."


def build_stress_results(stress_results_dir: Path, scenarios: List[str]) -> pd.DataFrame:
    det_rows = []
    abd_rows = []

    for scenario in scenarios:
        scenario_dir = stress_results_dir / scenario
        det_path = scenario_dir / "detection_metrics.tsv"
        abd_path = scenario_dir / "abundance_metrics.tsv"

        if not det_path.exists():
            print(f"WARNING: missing detection metrics for {scenario}: {det_path}", file=sys.stderr)
            continue

        det_rows.append(summarise_detection(det_path, scenario))
        abd_rows.append(summarise_abundance(abd_path, scenario))

    if not det_rows:
        raise RuntimeError(f"No detection metric files found under {stress_results_dir}")

    det = pd.concat(det_rows, ignore_index=True)
    abd = pd.concat(abd_rows, ignore_index=True) if abd_rows else pd.DataFrame()

    result = det.merge(abd, on=["test_id", "community_type"], how="left")

    baseline = result[result["test_id"] == "ST00_baseline"][
        ["community_type", "precision", "recall", "f1", "accuracy", "rmse", "bray_curtis"]
    ].copy()
    baseline = baseline.rename(
        columns={
            "precision": "baseline_precision",
            "recall": "baseline_recall",
            "f1": "baseline_f1",
            "accuracy": "baseline_accuracy",
            "rmse": "baseline_rmse",
            "bray_curtis": "baseline_bray_curtis",
        }
    )

    result = result.merge(baseline, on="community_type", how="left")
    result["delta_precision_vs_baseline"] = result["precision"] - result["baseline_precision"]
    result["delta_recall_vs_baseline"] = result["recall"] - result["baseline_recall"]
    result["delta_f1_vs_baseline"] = result["f1"] - result["baseline_f1"]
    result["delta_accuracy_vs_baseline"] = result["accuracy"] - result["baseline_accuracy"]
    result["delta_rmse_vs_baseline"] = result["rmse"] - result["baseline_rmse"]
    result["delta_bray_curtis_vs_baseline"] = result["bray_curtis"] - result["baseline_bray_curtis"]
    result["interpretation"] = result.apply(interpret_stress_row, axis=1)

    ordered = [
        "test_id",
        "description",
        "community_type",
        "precision",
        "recall",
        "f1",
        "accuracy",
        "rmse",
        "bray_curtis",
        "baseline_precision",
        "baseline_recall",
        "baseline_f1",
        "delta_precision_vs_baseline",
        "delta_recall_vs_baseline",
        "delta_f1_vs_baseline",
        "delta_rmse_vs_baseline",
        "delta_bray_curtis_vs_baseline",
        "interpretation",
    ]
    existing = [c for c in ordered if c in result.columns]
    return result[existing].sort_values(["test_id", "community_type"]).reset_index(drop=True)


def evaluate_ST03_invalid_counts(base_dir: Path) -> Dict[str, object]:
    path = base_dir / "technical_checks" / "ST03_invalid_count_table" / "asv_table.tsv"
    if not path.exists():
        return {
            "check_id": "ST03_invalid_count_table",
            "check_type": "invalid_count_table",
            "expected_behavior": "Non-numeric count values should be detected.",
            "observed_behavior": f"Missing file: {path}",
            "metric": "NA",
            "pass_fail": "FAIL",
            "notes": "Technical check file not found.",
        }

    df = pd.read_csv(path, sep="\t", dtype=str)
    sample_cols = list(df.columns[1:])
    numeric = df[sample_cols].apply(pd.to_numeric, errors="coerce")
    non_numeric_count = int(numeric.isna().sum().sum())

    return {
        "check_id": "ST03_invalid_count_table",
        "check_type": "invalid_count_table",
        "expected_behavior": "Non-numeric count values should be detected.",
        "observed_behavior": f"{non_numeric_count} non-numeric cells detected.",
        "metric": non_numeric_count,
        "pass_fail": "PASS" if non_numeric_count > 0 else "FAIL",
        "notes": "This verifies corrupt count tables can be detected before metric interpretation.",
    }


def evaluate_ST05_metadata_mismatch(base_dir: Path) -> Dict[str, object]:
    check_dir = base_dir / "technical_checks" / "ST05_metadata_mismatch"
    asv_path = check_dir / "asv_table.tsv"
    metadata_path = check_dir / "metadata_wrong.tsv"

    if not asv_path.exists() or not metadata_path.exists():
        return {
            "check_id": "ST05_metadata_mismatch",
            "check_type": "metadata_mismatch",
            "expected_behavior": "Metadata sample IDs should not match ASV table sample IDs.",
            "observed_behavior": f"Missing file(s): {asv_path} or {metadata_path}",
            "metric": "NA",
            "pass_fail": "FAIL",
            "notes": "Technical check file not found.",
        }

    asv = pd.read_csv(asv_path, sep="\t", nrows=1)
    asv_samples = set(map(str, asv.columns[1:]))

    meta = pd.read_csv(metadata_path, sep="\t", dtype=str)
    meta_col = find_col(meta, ["sample_id", "sample-id", "sample"]) or meta.columns[0]
    meta_samples = set(meta[meta_col].astype(str))

    matched = len(asv_samples & meta_samples)
    total_meta = len(meta_samples)
    total_asv = len(asv_samples)

    return {
        "check_id": "ST05_metadata_mismatch",
        "check_type": "metadata_mismatch",
        "expected_behavior": "Metadata sample IDs should not match ASV table sample IDs.",
        "observed_behavior": f"{matched} matching IDs out of {total_asv} ASV samples and {total_meta} metadata samples.",
        "metric": matched,
        "pass_fail": "PASS" if matched == 0 else "FAIL",
        "notes": "Expected 0 matches because metadata IDs contain an intentional suffix.",
    }


def evaluate_ST06a_contamination_blanks(base_dir: Path, min_blank_relative_abundance: float = 0.10) -> Dict[str, object]:
    scenario_dir = base_dir / "stress_inputs" / "ST06a_contamination_blanks"
    asv_path = scenario_dir / "asv_table.tsv"
    spike_path = scenario_dir / "spiked_contaminant_ASVs.tsv"
    metadata_path = base_dir / "metadata_PRJEB10949.tsv"

    if not asv_path.exists() or not spike_path.exists() or not metadata_path.exists():
        return {
            "check_id": "ST06a_contamination_blanks",
            "check_type": "blank_contamination_spike",
            "expected_behavior": "Spiked contaminant should be enriched in blank controls.",
            "observed_behavior": f"Missing file(s): {asv_path}, {spike_path}, or {metadata_path}",
            "metric": "NA",
            "pass_fail": "FAIL",
            "notes": "Contamination check file not found.",
        }

    asv = pd.read_csv(asv_path, sep="\t", index_col=0)
    spike = pd.read_csv(spike_path, sep="\t")
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)

    contam_col = find_col(spike, ["ASV_ID", "asv_id", "#OTU ID", "otu_id"]) or spike.columns[0]
    contam_asv = str(spike.iloc[0][contam_col])

    sample_col = find_col(metadata, ["sample_id", "sample-id", "sample"]) or metadata.columns[0]
    neg_col = find_col(metadata, ["is_negative_control", "negative_control"])
    community_col = find_col(metadata, ["community_type", "community"])

    if neg_col:
        blank_samples = metadata[
            metadata[neg_col].astype(str).str.lower().isin(["yes", "true", "1"])
        ][sample_col].tolist()
    elif community_col:
        blank_samples = metadata[
            metadata[community_col].astype(str).str.lower().str.contains("negative|blank|h2o", regex=True)
        ][sample_col].tolist()
    else:
        blank_samples = []

    blank_samples = [s for s in blank_samples if s in asv.columns]

    if contam_asv not in asv.index:
        return {
            "check_id": "ST06a_contamination_blanks",
            "check_type": "blank_contamination_spike",
            "expected_behavior": "Spiked contaminant should be present in ASV table.",
            "observed_behavior": f"Contaminant ASV {contam_asv} not found in ASV table.",
            "metric": "NA",
            "pass_fail": "FAIL",
            "notes": "Spike metadata and ASV table are inconsistent.",
        }

    if not blank_samples:
        return {
            "check_id": "ST06a_contamination_blanks",
            "check_type": "blank_contamination_spike",
            "expected_behavior": "Blank samples should be identifiable from metadata.",
            "observed_behavior": "No blank samples found in metadata.",
            "metric": "NA",
            "pass_fail": "FAIL",
            "notes": "Check metadata_PRJEB10949.tsv.",
        }

    totals = asv[blank_samples].sum(axis=0)
    rel = asv.loc[contam_asv, blank_samples] / totals.replace(0, pd.NA)
    mean_rel = float(rel.mean(skipna=True))

    return {
        "check_id": "ST06a_contamination_blanks",
        "check_type": "blank_contamination_spike",
        "expected_behavior": f"Spiked contaminant relative abundance in blanks should be >= {min_blank_relative_abundance:.0%}.",
        "observed_behavior": f"Mean relative abundance of {contam_asv} in blanks = {mean_rel:.3f}.",
        "metric": mean_rel,
        "pass_fail": "PASS" if mean_rel >= min_blank_relative_abundance else "FAIL",
        "notes": "This scenario is intended for contamination filtering/decontam rather than validation_metrics.py.",
    }


def build_technical_results(base_dir: Path) -> pd.DataFrame:
    rows = [
        evaluate_ST03_invalid_counts(base_dir),
        evaluate_ST05_metadata_mismatch(base_dir),
        evaluate_ST06a_contamination_blanks(base_dir),
    ]
    return pd.DataFrame(rows)


def fmt_value(v: object) -> str:
    try:
        if pd.isna(v):
            return "NA"
        if isinstance(v, float):
            return f"{v:.4f}"
        return str(v)
    except Exception:
        return str(v)


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows available._"

    cols = list(df.columns)
    lines = []
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(fmt_value(row[c]) for c in cols) + " |")
    return "\n".join(lines)


def write_report(stress_results: pd.DataFrame, technical_results: pd.DataFrame, outpath: Path) -> None:
    main_cols = [
        "test_id",
        "community_type",
        "precision",
        "recall",
        "f1",
        "rmse",
        "bray_curtis",
        "delta_f1_vs_baseline",
        "interpretation",
    ]
    main_cols = [c for c in main_cols if c in stress_results.columns]

    technical_cols = ["check_id", "check_type", "observed_behavior", "pass_fail", "notes"]
    technical_cols = [c for c in technical_cols if c in technical_results.columns]

    lines = [
        "# PRJEB10949 Stress Test Report",
        "",
        "## Objective",
        "",
        "This report summarises automated stress testing of the Group C validation workflow using the PRJEB10949 BEI mock-community dataset. The main scenarios are compared against the ST00 baseline using detection and abundance metrics generated by `validation_metrics.py`.",
        "",
        "## Main quantitative stress-test results",
        "",
        dataframe_to_markdown(stress_results[main_cols]),
        "",
        "## Technical and contamination checks",
        "",
        dataframe_to_markdown(technical_results[technical_cols]),
        "",
        "## Interpretation notes",
        "",
        "- `delta_f1_vs_baseline` is calculated within each community type using ST00 as reference.",
        "- ST03 and ST05 are technical checks and are not expected to produce F1/recall metrics.",
        "- ST06a targets blank controls and is intended for contamination filtering evaluation.",
        "- ST06b targets biological mock samples and is included in the quantitative validation summary.",
        "",
        "## Generated files",
        "",
        "- `stress_test_results.tsv`",
        "- `technical_check_results.tsv`",
        "- `stress_test_report.md`",
        "",
    ]
    outpath.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Summarise PRJEB10949 stress-test outputs.")
    ap.add_argument("--base-dir", default=".", help="Path to stress_tests directory. Default: current directory.")
    ap.add_argument("--stress-results-dir", default="stress_results", help="Directory containing validation_metrics.py outputs.")
    ap.add_argument("--out-stress", default="stress_test_results.tsv", help="Output TSV for quantitative stress-test summary.")
    ap.add_argument("--out-technical", default="technical_check_results.tsv", help="Output TSV for technical/contamination check summary.")
    ap.add_argument("--out-report", default="stress_test_report.md", help="Output Markdown report.")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    base_dir = Path(args.base_dir).resolve()
    stress_results_dir = base_dir / args.stress_results_dir

    if not stress_results_dir.exists():
        raise FileNotFoundError(
            f"Stress results directory not found: {stress_results_dir}\n"
            "Run `bash scripts/run_validation_metrics_all.sh` first."
        )

    stress_results = build_stress_results(stress_results_dir, MAIN_SCENARIOS)
    technical_results = build_technical_results(base_dir)

    stress_out = base_dir / args.out_stress
    technical_out = base_dir / args.out_technical
    report_out = base_dir / args.out_report

    stress_results.to_csv(stress_out, sep="\t", index=False)
    technical_results.to_csv(technical_out, sep="\t", index=False)
    write_report(stress_results, technical_results, report_out)

    print(f"✓ Wrote {stress_out}")
    print(f"✓ Wrote {technical_out}")
    print(f"✓ Wrote {report_out}")


if __name__ == "__main__":
    main()
