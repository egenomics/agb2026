"""
parsers.py — QIIME2 TSV parsing functions for the MicroSee report generator.

Pure parsing functions — no FastAPI, no HTTP, no I/O.
Each function takes a raw string and returns a typed result or raises ValueError.

The integration step (joining all parsed results) lives in integrator.py.
It is re-exported here so callers can use either import path.
"""

from __future__ import annotations

import io
import logging
import re

import numpy as np
import pandas as pd

from .integrator import integrate
from .models import (
    AlphaDiversityEntry,
    AlphaDiversityResult,
    DistanceMatrixResult,
    FeatureTableResult,
    MetadataResult,
    SampleMetadata,
    TaxonomyResult,
)

logger = logging.getLogger(__name__)

__all__ = [
    "integrate",
    "parse_alpha_diversity",
    "parse_distance_matrix",
    "parse_feature_table",
    "parse_metadata",
    "parse_taxonomy",
]


# ── Shared helpers ────────────────────────────────────────────────────────────


def _read_tsv(content: str, skip_prefixes: tuple[str, ...] = ("#",)) -> pd.DataFrame:
    """
    Parse TSV content into a DataFrame, stripping QIIME2 comment lines.
    Raises ValueError if the content is empty or cannot be parsed.
    """
    lines = [
        line
        for line in content.strip().splitlines()
        if line.strip() and not any(line.startswith(p) for p in skip_prefixes)
    ]
    if not lines:
        raise ValueError("File is empty or contains only comment lines.")

    try:
        return pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", low_memory=False)
    except Exception as exc:
        raise ValueError(f"Could not parse TSV: {exc}") from exc


def _find_column(df: pd.DataFrame, patterns: list[str]) -> str | None:
    """Return the first column matching any of the regex patterns (case-insensitive)."""
    for pattern in patterns:
        for col in df.columns:
            if re.search(pattern, col, re.IGNORECASE):
                return col
    return None


def _require_column(df: pd.DataFrame, patterns: list[str], label: str) -> str:
    """Like _find_column but raises ValueError if nothing matches."""
    col = _find_column(df, patterns)
    if col is None:
        raise ValueError(
            f"Could not find a '{label}' column. "
            f"Expected a column matching: {patterns}. "
            f"Found columns: {list(df.columns)}"
        )
    return col


# ── Feature table ─────────────────────────────────────────────────────────────


def parse_feature_table(content: str) -> FeatureTableResult:
    """
    Parse a QIIME2 feature-table.tsv export.

    Expected format:
        # Constructed from biom file
        #OTU ID    sample1    sample2    ...
        asv001     100        0          ...

    Returns FeatureTableResult with samples, features, and count matrix.
    Raises ValueError on format errors.
    """
    lines = [
        line
        for line in content.strip().splitlines()
        if line.strip() and not line.startswith("# Constructed")
    ]
    if not lines:
        raise ValueError("Feature table is empty.")

    content_clean = "\n".join(lines)
    try:
        df = pd.read_csv(io.StringIO(content_clean), sep="\t", index_col=0)
    except Exception as exc:
        raise ValueError(f"Cannot parse feature table: {exc}") from exc

    df = df.dropna(how="all").fillna(0)

    if df.empty:
        raise ValueError("Feature table has no data rows.")
    if len(df.columns) == 0:
        raise ValueError("Feature table has no sample columns.")

    try:
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    except Exception as exc:
        raise ValueError("Feature table contains non-numeric count values.") from exc

    features: list[str] = [str(f) for f in df.index]
    samples: list[str] = [str(s) for s in df.columns]
    counts: dict[str, dict[str, float]] = {
        feat: {samp: float(str(df.at[feat, samp])) for samp in samples} for feat in features
    }

    logger.info("Parsed feature table: %d features × %d samples", len(features), len(samples))

    return FeatureTableResult(
        samples=samples,
        features=features,
        counts=counts,
        n_samples=len(samples),
        n_features=len(features),
    )


# ── Taxonomy ──────────────────────────────────────────────────────────────────

_FAMILY_RE = re.compile(r"f__([^;]+)")
_SKIP_VALUES = {"", "uncultured", "unknown", "unidentified", "metagenome"}


def _extract_family(taxon: object) -> str:
    """Extract family-level name from a QIIME2/SILVA taxonomy string."""
    if not isinstance(taxon, str):
        return "Unclassified"
    m = _FAMILY_RE.search(taxon)
    if m:
        name = m.group(1).strip()
        if name.lower() not in _SKIP_VALUES:
            return name
    return "Unclassified"


def parse_taxonomy(content: str) -> TaxonomyResult:
    """
    Parse a QIIME2 taxonomy.tsv export.

    Expected format:
        Feature ID    Taxon    Confidence
        asv001        d__Bact;...;f__Bacteroidaceae;...    0.99

    Returns TaxonomyResult mapping feature_id → family name.
    Raises ValueError on format errors.
    """
    df = _read_tsv(content)

    fid_col = _require_column(df, [r"^feature.?id$", r"^id$", r"^#?otu.?id$"], "Feature ID")
    tax_col = _require_column(df, [r"taxon", r"taxonomy", r"classif"], "Taxon")

    assignments: dict[str, str] = {}
    for _, row in df.iterrows():
        fid = str(row[fid_col]).strip()
        family = _extract_family(str(row.get(tax_col, "")))
        assignments[fid] = family

    if not assignments:
        raise ValueError("Taxonomy file has no valid rows.")

    unclassified = sum(1 for v in assignments.values() if v == "Unclassified")
    unclassified_pct = round(unclassified / len(assignments) * 100, 1)
    unique_families = sorted(set(v for v in assignments.values() if v != "Unclassified"))

    if unclassified_pct > 80:
        logger.warning(
            "%.1f%% of features are unclassified at family level. "
            "Check that your taxonomy file uses SILVA or Greengenes format.",
            unclassified_pct,
        )

    logger.info(
        "Parsed taxonomy: %d features, %d families, %.1f%% unclassified",
        len(assignments),
        len(unique_families),
        unclassified_pct,
    )

    return TaxonomyResult(
        assignments=assignments,
        unclassified_pct=unclassified_pct,
        unique_families=unique_families,
    )


# ── Metadata ──────────────────────────────────────────────────────────────────


def _parse_time_days(timepoint_str: str) -> int | None:
    """
    Convert a timepoint string to integer days.
    Examples: 'T0'→0, 'T84'→84, 'Week12'→84, 'baseline'→0, '0'→0
    Returns None if unparseable.
    """
    s = str(timepoint_str).strip().lower()
    if re.fullmatch(r"\d+", s):
        return int(s)
    m = re.match(r"t(\d+)", s)
    if m:
        return int(m.group(1))
    m = re.match(r"(?:week|wk)(\d+)", s)
    if m:
        return int(m.group(1)) * 7
    if s in ("baseline", "pre", "week0", "wk0", "visit0"):
        return 0
    return None


def parse_metadata(content: str) -> MetadataResult:
    """
    Parse a QIIME2 metadata.tsv export.

    Required column: sample-id (or sampleid, #SampleID, etc.)
    Recognised columns (case-insensitive):
        group / treatment / condition / intervention
        timepoint / time_point / visit / week
        patient / subject / individual / participant
        sixmwt / 6mwt / walk_distance
        il18 / il-18 / interleukin18

    Returns MetadataResult.
    Raises ValueError if sample-id column is missing.
    """
    lines = [
        line
        for line in content.strip().splitlines()
        if line.strip() and not line.startswith("#q2:types")
    ]
    df = pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", dtype=str).fillna("")

    sid_col = _require_column(df, [r"^#?sample[-_]?id$", r"^samplename$", r"^id$"], "sample-id")
    grp_col = _find_column(df, [r"^group$", r"^treatment$", r"^condition$", r"^intervention$"])
    tp_col = _find_column(df, [r"timepoint", r"time[-_]?point", r"^visit$", r"^week$"])
    pat_col = _find_column(df, [r"^patient$", r"^subject$", r"^individual$", r"^participant$"])
    mwt_col = _find_column(df, [r"6mwt", r"sixmwt", r"walk"])
    il18_col = _find_column(df, [r"il.?18", r"interleukin"])

    samples: list[SampleMetadata] = []

    for _, row in df.iterrows():
        sid = str(row[sid_col]).strip()
        if not sid:
            continue

        grp_raw = str(row[grp_col]).strip() if grp_col else ""
        tp_raw = str(row[tp_col]).strip() if tp_col else ""
        pat_raw = str(row[pat_col]).strip() if pat_col else re.sub(r"_T\d+$", "", sid)

        base_group = grp_raw or sid
        time_days = _parse_time_days(tp_raw) if tp_raw else None
        full_group = f"{base_group}_{tp_raw}" if (base_group and tp_raw) else base_group

        sixmwt = 0.0
        il18 = 0.0
        if mwt_col:
            try:
                sixmwt = float(str(row[mwt_col]))
            except (ValueError, TypeError):
                pass
        if il18_col:
            try:
                il18 = float(str(row[il18_col]))
            except (ValueError, TypeError):
                pass

        samples.append(
            SampleMetadata(
                sample_id=sid,
                group=full_group,
                base_group=base_group,
                timepoint=tp_raw,
                time=time_days,
                patient=pat_raw,
                sixmwt=sixmwt,
                il18=il18,
            )
        )

    if not samples:
        raise ValueError("Metadata file has no valid sample rows.")

    has_clinical = any(s.sixmwt > 0 or s.il18 > 0 for s in samples)
    groups = sorted(set(s.group for s in samples))
    base_groups = sorted(set(s.base_group for s in samples))
    timepoints = sorted(set(s.timepoint for s in samples if s.timepoint))

    logger.info(
        "Parsed metadata: %d samples, %d groups, clinical=%s",
        len(samples),
        len(groups),
        "yes" if has_clinical else "no",
    )

    return MetadataResult(
        samples=samples,
        groups=groups,
        base_groups=base_groups,
        timepoints=timepoints,
        has_clinical=has_clinical,
        n_samples=len(samples),
    )


# ── Alpha diversity ───────────────────────────────────────────────────────────

_ALPHA_PATTERNS: dict[str, list[str]] = {
    "shannon": [r"shannon"],
    "simpson": [r"simpson"],
    "observed": [r"observed"],
    "faith_pd": [r"faith"],
    "pielou": [r"pielou", r"evenness"],
}


def _safe_float(row: pd.Series, col: str | None) -> float:
    if col is None or col not in row:
        return 0.0
    try:
        return round(float(str(row[col])), 4)
    except (ValueError, TypeError):
        return 0.0


def parse_alpha_diversity(content: str) -> AlphaDiversityResult:
    """
    Parse a QIIME2 alpha-diversity.tsv export.
    Handles multiple metrics in one file (QIIME2 exports one metric per file,
    but we accept merged files too).

    Returns AlphaDiversityResult.
    Raises ValueError if sample-id column is missing.
    """
    df = _read_tsv(content)
    df.columns = [c.strip() for c in df.columns]

    sid_col = _require_column(df, [r"^#?sample[-_]?id$", r"^id$"], "sample-id")

    detected: dict[str, str] = {}
    for metric, patterns in _ALPHA_PATTERNS.items():
        col = _find_column(df, patterns)
        if col and col != sid_col:
            detected[metric] = col

    if not detected:
        raise ValueError(
            f"No recognised alpha diversity metrics found. "
            f"Expected columns like: shannon_entropy, observed_features, faith_pd. "
            f"Found: {list(df.columns)}"
        )

    entries: list[AlphaDiversityEntry] = []
    for _, row in df.iterrows():
        sid = str(row[sid_col]).strip()
        if not sid:
            continue

        entries.append(
            AlphaDiversityEntry(
                sample_id=sid,
                shannon=_safe_float(row, detected.get("shannon")),
                simpson=_safe_float(row, detected.get("simpson")),
                observed=_safe_float(row, detected.get("observed")),
                faith_pd=_safe_float(row, detected.get("faith_pd")),
                pielou=_safe_float(row, detected.get("pielou")),
            )
        )

    if not entries:
        raise ValueError("Alpha diversity file has no valid sample rows.")

    logger.info(
        "Parsed alpha diversity: %d samples, metrics: %s",
        len(entries),
        list(detected.keys()),
    )

    return AlphaDiversityResult(
        samples=entries,
        metrics_present=list(detected.keys()),
    )


# ── Distance matrix ───────────────────────────────────────────────────────────


def parse_distance_matrix(content: str) -> DistanceMatrixResult:
    """
    Parse a QIIME2 distance-matrix.tsv export.

    Format (note leading tab on header line — QIIME2 convention):
        \\t    sample1    sample2    ...
        sample1    0       0.42       ...
        sample2    0.42    0          ...

    Returns DistanceMatrixResult.
    Raises ValueError on format or symmetry errors.
    """
    lines = [line for line in content.splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError("Distance matrix file has fewer than 2 rows.")

    try:
        df = pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", index_col=0)
    except Exception as exc:
        raise ValueError(f"Cannot parse distance matrix: {exc}") from exc

    if df.empty:
        raise ValueError("Distance matrix is empty after parsing.")

    samples = df.index.astype(str).tolist()
    n = len(samples)

    if df.shape[0] != df.shape[1]:
        raise ValueError(
            f"Distance matrix is not square: {df.shape[0]} rows × {df.shape[1]} columns."
        )

    try:
        mat = df.astype(float).values
    except Exception as exc:
        raise ValueError(f"Distance matrix contains non-numeric values: {exc}") from exc

    max_asymmetry = float(np.max(np.abs(mat - mat.T)))
    if max_asymmetry > 1e-4:
        logger.warning("Distance matrix asymmetry detected: max=%.6f", max_asymmetry)

    logger.info("Parsed distance matrix: %d × %d", n, n)

    return DistanceMatrixResult(samples=samples, matrix=mat.tolist(), n=n)
