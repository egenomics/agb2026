"""
parsers.py — QIIME2 TSV parsing functions for the MicroSee report generator.

Pure parsing functions — no FastAPI, no HTTP, no I/O.
Each function takes a raw string and returns a typed result or raises ValueError.
"""

from __future__ import annotations

import re
import io
import logging
from typing import Optional

import pandas as pd
import numpy as np

from .models import (
    FeatureTableResult,
    TaxonomyResult,
    MetadataResult,
    SampleMetadata,
    AlphaDiversityResult,
    AlphaDiversityEntry,
    DistanceMatrixResult,
    IntegrateResult,
    SampleRow,
)

logger = logging.getLogger(__name__)


# ── Shared helpers ────────────────────────────────────────────────────────────

def _read_tsv(content: str, skip_prefixes: tuple[str, ...] = ("#",)) -> pd.DataFrame:
    """
    Parse TSV content into a DataFrame, stripping QIIME2 comment lines.
    Raises ValueError if the content is empty or cannot be parsed.
    """
    lines = [
        line for line in content.strip().splitlines()
        if line.strip() and not any(line.startswith(p) for p in skip_prefixes)
    ]
    if not lines:
        raise ValueError("File is empty or contains only comment lines.")

    try:
        return pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", low_memory=False)
    except Exception as exc:
        raise ValueError(f"Could not parse TSV: {exc}") from exc


def _find_column(df: pd.DataFrame, patterns: list[str]) -> Optional[str]:
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
    # Strip the biom comment line specifically
    lines = [
        line for line in content.strip().splitlines()
        if line.strip() and not line.startswith("# Constructed")
    ]
    if not lines:
        raise ValueError("Feature table is empty.")

    # QIIME2 sometimes uses '#OTU ID' as the index column name
    content_clean = "\n".join(lines)
    try:
        df = pd.read_csv(io.StringIO(content_clean), sep="\t", index_col=0)
    except Exception as exc:
        raise ValueError(f"Cannot parse feature table: {exc}") from exc

    # Drop empty rows/columns
    df = df.dropna(how="all").fillna(0)

    if df.empty:
        raise ValueError("Feature table has no data rows.")
    if len(df.columns) == 0:
        raise ValueError("Feature table has no sample columns.")

    # Ensure all values are numeric
    try:
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    except Exception:
        raise ValueError("Feature table contains non-numeric count values.")

    features = df.index.astype(str).tolist()
    samples  = df.columns.astype(str).tolist()
    counts   = {
        feat: {samp: float(df.at[feat, samp]) for samp in samples}
        for feat in features
    }

    logger.info(f"Parsed feature table: {len(features)} features × {len(samples)} samples")

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


def _extract_family(taxon: str) -> str:
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
        fid    = str(row[fid_col]).strip()
        family = _extract_family(str(row.get(tax_col, "")))
        assignments[fid] = family

    if not assignments:
        raise ValueError("Taxonomy file has no valid rows.")

    unclassified     = sum(1 for v in assignments.values() if v == "Unclassified")
    unclassified_pct = round(unclassified / len(assignments) * 100, 1)
    unique_families  = sorted(set(v for v in assignments.values() if v != "Unclassified"))

    if unclassified_pct > 80:
        logger.warning(
            f"{unclassified_pct}% of features are unclassified at family level. "
            "Check that your taxonomy file uses SILVA or Greengenes format."
        )

    logger.info(
        f"Parsed taxonomy: {len(assignments)} features, "
        f"{len(unique_families)} families, "
        f"{unclassified_pct}% unclassified"
    )

    return TaxonomyResult(
        assignments=assignments,
        unclassified_pct=unclassified_pct,
        unique_families=unique_families,
    )


# ── Metadata ──────────────────────────────────────────────────────────────────

def _parse_time_days(timepoint_str: str) -> Optional[int]:
    """
    Convert a timepoint string to integer days.
    Examples: 'T0'→0, 'T84'→84, 'Week12'→84, 'baseline'→0, '0'→0
    Returns None if unparseable.
    """
    s = str(timepoint_str).strip().lower()
    # Direct numeric
    if re.fullmatch(r"\d+", s):
        return int(s)
    # T<N> pattern
    m = re.match(r"t(\d+)", s)
    if m:
        return int(m.group(1))
    # Week<N> or wk<N> → convert to days
    m = re.match(r"(?:week|wk)(\d+)", s)
    if m:
        return int(m.group(1)) * 7
    # Baseline / control synonyms
    if s in ("baseline", "t0", "pre", "week0", "wk0", "visit0"):
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
    # Skip QIIME2 type directive line (#q2:types)
    lines = [
        l for l in content.strip().splitlines()
        if l.strip() and not l.startswith("#q2:types")
    ]
    df = pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", dtype=str).fillna("")

    sid_col = _require_column(
        df,
        [r"^#?sample[-_]?id$", r"^samplename$", r"^id$"],
        "sample-id"
    )
    grp_col  = _find_column(df, [r"^group$", r"^treatment$", r"^condition$", r"^intervention$"])
    tp_col   = _find_column(df, [r"timepoint", r"time[-_]?point", r"^visit$", r"^week$"])
    pat_col  = _find_column(df, [r"^patient$", r"^subject$", r"^individual$", r"^participant$"])
    mwt_col  = _find_column(df, [r"6mwt", r"sixmwt", r"walk"])
    il18_col = _find_column(df, [r"il.?18", r"interleukin"])

    samples: list[SampleMetadata] = []

    for _, row in df.iterrows():
        sid = str(row[sid_col]).strip()
        if not sid:
            continue

        grp_raw = str(row[grp_col]).strip()  if grp_col  else ""
        tp_raw  = str(row[tp_col]).strip()   if tp_col   else ""
        pat_raw = str(row[pat_col]).strip()  if pat_col  else re.sub(r"_T\d+$", "", sid)

        # Derive base_group and full group
        base_group = grp_raw or sid
        time_days  = _parse_time_days(tp_raw) if tp_raw else None
        full_group = f"{base_group}_{tp_raw}" if (base_group and tp_raw) else base_group

        # Clinical columns
        sixmwt = 0.0
        il18   = 0.0
        if mwt_col:
            try:   sixmwt = float(row[mwt_col])
            except (ValueError, TypeError): pass
        if il18_col:
            try:   il18 = float(row[il18_col])
            except (ValueError, TypeError): pass

        samples.append(SampleMetadata(
            sample_id=sid,
            group=full_group,
            base_group=base_group,
            timepoint=tp_raw,
            time=time_days,
            patient=pat_raw,
            sixmwt=sixmwt,
            il18=il18,
        ))

    if not samples:
        raise ValueError("Metadata file has no valid sample rows.")

    has_clinical  = any(s.sixmwt > 0 or s.il18 > 0 for s in samples)
    groups        = sorted(set(s.group      for s in samples))
    base_groups   = sorted(set(s.base_group for s in samples))
    timepoints    = sorted(set(s.timepoint  for s in samples if s.timepoint))

    logger.info(
        f"Parsed metadata: {len(samples)} samples, "
        f"{len(groups)} groups, clinical={'yes' if has_clinical else 'no'}"
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
    "shannon":  [r"shannon"],
    "simpson":  [r"simpson"],
    "observed": [r"observed"],
    "faith_pd": [r"faith"],
    "pielou":   [r"pielou", r"evenness"],
}


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

    # Detect which metric columns exist
    detected: dict[str, str] = {}  # metric_name → actual column name
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

        def _safe_float(col: Optional[str]) -> float:
            if col is None or col not in row:
                return 0.0
            try:   return round(float(row[col]), 4)
            except: return 0.0

        entries.append(AlphaDiversityEntry(
            sample_id=sid,
            shannon  =_safe_float(detected.get("shannon")),
            simpson  =_safe_float(detected.get("simpson")),
            observed =_safe_float(detected.get("observed")),
            faith_pd =_safe_float(detected.get("faith_pd")),
            pielou   =_safe_float(detected.get("pielou")),
        ))

    if not entries:
        raise ValueError("Alpha diversity file has no valid sample rows.")

    logger.info(
        f"Parsed alpha diversity: {len(entries)} samples, "
        f"metrics: {list(detected.keys())}"
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

    Uses pandas with index_col=0 to correctly handle the leading tab.
    Returns DistanceMatrixResult.
    Raises ValueError on format or symmetry errors.
    """
    lines = [l for l in content.splitlines() if l.strip()]
    if len(lines) < 2:
        raise ValueError("Distance matrix file has fewer than 2 rows.")

    try:
        df = pd.read_csv(
            io.StringIO("\n".join(lines)),
            sep="\t",
            index_col=0,
        )
    except Exception as exc:
        raise ValueError(f"Cannot parse distance matrix: {exc}") from exc

    if df.empty:
        raise ValueError("Distance matrix is empty after parsing.")

    samples = df.index.astype(str).tolist()
    n       = len(samples)

    if df.shape[0] != df.shape[1]:
        raise ValueError(
            f"Distance matrix is not square: "
            f"{df.shape[0]} rows × {df.shape[1]} columns."
        )

    try:
        mat = df.astype(float).values
    except Exception as exc:
        raise ValueError(f"Distance matrix contains non-numeric values: {exc}") from exc

    # Symmetry check
    max_asymmetry = float(np.max(np.abs(mat - mat.T)))
    if max_asymmetry > 1e-4:
        logger.warning(f"Distance matrix asymmetry detected: max={max_asymmetry:.6f}")

    logger.info(f"Parsed distance matrix: {n} × {n}")

    return DistanceMatrixResult(
        samples=samples,
        matrix=mat.tolist(),
        n=n,
    )


# ── Integration ───────────────────────────────────────────────────────────────

def integrate(
    feature_table: FeatureTableResult,
    taxonomy:      TaxonomyResult,
    metadata:      MetadataResult,
    alpha:         Optional[AlphaDiversityResult] = None,
) -> IntegrateResult:
    """
    Combine parsed QIIME2 outputs into chart-ready SampleRow objects.

    Steps:
    1. Aggregate ASV counts to family level using taxonomy assignments
    2. Convert to relative abundance (%)
    3. Join with metadata
    4. Join with alpha diversity (if provided)
    5. Compute Shannon + Simpson if not in alpha file
    6. Sort taxa by mean abundance descending

    Returns IntegrateResult.
    Raises ValueError if no samples survive the join.
    """
    warnings: list[str] = []

    # ── Step 1: aggregate to family level ──
    family_counts: dict[str, dict[str, float]] = {}  # family → {sample_id: count}

    for feature_id, sample_counts in feature_table.counts.items():
        family = taxonomy.assignments.get(feature_id, "Unclassified")
        if family == "Unclassified":
            continue
        if family not in family_counts:
            family_counts[family] = {s: 0.0 for s in feature_table.samples}
        for sample_id, count in sample_counts.items():
            family_counts[family][sample_id] = family_counts[family].get(sample_id, 0.0) + count

    if not family_counts:
        raise ValueError(
            "No features could be assigned to a family. "
            "Check that taxonomy.tsv uses the same feature IDs as feature-table.tsv "
            "and uses SILVA/Greengenes format with 'f__' family prefixes."
        )

    # ── Step 2: relative abundance per sample ──
    families = list(family_counts.keys())
    rel_ab: dict[str, dict[str, float]] = {}  # sample_id → {family: pct}

    for sample_id in feature_table.samples:
        totals = {fam: family_counts[fam].get(sample_id, 0.0) for fam in families}
        total  = sum(totals.values()) or 1.0
        rel_ab[sample_id] = {fam: round(cnt / total * 100, 3) for fam, cnt in totals.items()}

    # ── Step 3: build alpha lookup ──
    alpha_lookup: dict[str, AlphaDiversityEntry] = {}
    if alpha:
        alpha_lookup = {e.sample_id: e for e in alpha.samples}

    # ── Step 4: build metadata lookup ──
    meta_lookup: dict[str, SampleMetadata] = {s.sample_id: s for s in metadata.samples}

    # ── Step 5: join everything ──
    rows: list[SampleRow] = []
    missing_meta: list[str] = []

    for sample_id in feature_table.samples:
        meta = meta_lookup.get(sample_id)
        if meta is None:
            missing_meta.append(sample_id)
            continue

        ab    = rel_ab.get(sample_id, {})
        alph  = alpha_lookup.get(sample_id)

        # Compute Shannon/Simpson if not in alpha file
        if alph and alph.shannon > 0:
            shannon = alph.shannon
            simpson = alph.simpson
        else:
            props   = np.array([ab.get(fam, 0.0) for fam in families], dtype=float)
            tot     = props.sum() or 1.0
            p       = props / tot
            shannon = float(-np.sum(p[p > 0] * np.log(p[p > 0])))
            simpson = float(1 - np.sum(p ** 2))

        row = SampleRow(
            sample_id  = sample_id,
            patient    = meta.patient,
            group      = meta.group,
            base_group = meta.base_group,
            timepoint  = meta.timepoint,
            time       = meta.time,
            shannon    = round(shannon, 4),
            simpson    = round(simpson, 4),
            pielou     = round(alph.pielou,   4) if alph else 0.0,
            observed   = round(alph.observed, 4) if alph else 0.0,
            faith_pd   = round(alph.faith_pd, 4) if alph else 0.0,
            sixmwt     = meta.sixmwt,
            il18       = meta.il18,
            **ab,  # taxon abundances as extra fields
        )
        rows.append(row)

    if not rows:
        raise ValueError(
            f"No samples survived the metadata join. "
            f"Feature table samples: {feature_table.samples[:5]}. "
            f"Metadata samples: {[s.sample_id for s in metadata.samples[:5]]}. "
            "Make sure sample IDs match exactly between files."
        )

    if missing_meta:
        warnings.append(
            f"{len(missing_meta)} samples in feature-table.tsv have no metadata: "
            f"{missing_meta[:5]}"
            f"{'...' if len(missing_meta) > 5 else ''}"
        )

    # ── Step 6: sort taxa by mean abundance descending ──
    mean_ab    = {fam: np.mean([rel_ab[s].get(fam, 0.0) for s in rel_ab]) for fam in families}
    taxa_sorted = sorted(families, key=lambda f: mean_ab[f], reverse=True)

    has_clinical = any(r.sixmwt > 0 or r.il18 > 0 for r in rows)
    groups       = sorted(set(r.group for r in rows))

    logger.info(
        f"Integration complete: {len(rows)} samples, "
        f"{len(taxa_sorted)} families, "
        f"groups={groups}"
    )

    return IntegrateResult(
        rows         = rows,
        taxa         = taxa_sorted,
        n_samples    = len(rows),
        n_taxa       = len(taxa_sorted),
        groups       = groups,
        has_clinical = has_clinical,
        warnings     = warnings,
    )
