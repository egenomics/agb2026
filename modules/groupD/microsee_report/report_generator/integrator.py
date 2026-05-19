"""
integrator.py — joins all parsed QIIME2 outputs into chart-ready SampleRow objects.

Kept separate from parsers.py because this step combines already-parsed data
rather than converting raw text, and contains the most complex business logic.
"""

from __future__ import annotations

import logging

import numpy as np

from .models import (
    AlphaDiversityEntry,
    AlphaDiversityResult,
    FeatureTableResult,
    IntegrateResult,
    MetadataResult,
    SampleRow,
    TaxonomyResult,
)

logger = logging.getLogger(__name__)


def integrate(
    feature_table: FeatureTableResult,
    taxonomy: TaxonomyResult,
    metadata: MetadataResult,
    alpha: AlphaDiversityResult | None = None,
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
    family_counts: dict[str, dict[str, float]] = {}

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
    rel_ab: dict[str, dict[str, float]] = {}

    for sample_id in feature_table.samples:
        totals = {fam: family_counts[fam].get(sample_id, 0.0) for fam in families}
        total = sum(totals.values()) or 1.0
        rel_ab[sample_id] = {fam: round(cnt / total * 100, 3) for fam, cnt in totals.items()}

    # ── Step 3: build lookups ──
    alpha_lookup: dict[str, AlphaDiversityEntry] = (
        {e.sample_id: e for e in alpha.samples} if alpha else {}
    )
    meta_lookup = {s.sample_id: s for s in metadata.samples}

    if alpha is None:
        logger.warning(
            "No alpha diversity file supplied — Shannon/Simpson will be estimated from "
            "family-level abundances (lower than ASV-level QIIME2 metrics). "
            "Pass --alpha with the merged alpha-diversity.tsv export."
        )

    # ── Step 4: join everything ──
    rows: list[SampleRow] = []
    missing_meta: list[str] = []

    for sample_id in feature_table.samples:
        meta = meta_lookup.get(sample_id)
        if meta is None:
            missing_meta.append(sample_id)
            continue

        ab = rel_ab.get(sample_id, {})
        alph = alpha_lookup.get(sample_id)

        # Compute Shannon/Simpson if not in alpha file.
        # WARNING: fallback computes from family-level relative abundances, which
        # collapses ASV-level variation and systematically underestimates diversity.
        # Pass --alpha with the QIIME2 alpha-diversity export to get accurate values.
        if alph and alph.shannon > 0:
            shannon = alph.shannon
            simpson = alph.simpson
        else:
            props = np.array([ab.get(fam, 0.0) for fam in families], dtype=float)
            tot = props.sum() or 1.0
            p = props / tot
            shannon = float(-np.sum(p[p > 0] * np.log(p[p > 0])))
            simpson = float(1 - np.sum(p**2))

        row = SampleRow(
            sample_id=sample_id,
            patient=meta.patient,
            group=meta.group,
            base_group=meta.base_group,
            timepoint=meta.timepoint,
            time=meta.time,
            shannon=round(shannon, 4),
            simpson=round(simpson, 4),
            pielou=round(alph.pielou, 4) if alph else 0.0,
            observed=round(alph.observed, 4) if alph else 0.0,
            faith_pd=round(alph.faith_pd, 4) if alph else 0.0,
            sixmwt=meta.sixmwt,
            il18=meta.il18,
            **ab,
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

    # ── Step 5: sort taxa by mean abundance descending ──
    mean_ab = {fam: np.mean([rel_ab[s].get(fam, 0.0) for s in rel_ab]) for fam in families}
    taxa_sorted = sorted(families, key=lambda f: mean_ab[f], reverse=True)

    has_clinical = any(r.sixmwt > 0 or r.il18 > 0 for r in rows)
    groups = sorted(set(r.group for r in rows))

    logger.info(
        "Integration complete: %d samples, %d families, groups=%s",
        len(rows),
        len(taxa_sorted),
        groups,
    )

    return IntegrateResult(
        rows=rows,
        taxa=taxa_sorted,
        n_samples=len(rows),
        n_taxa=len(taxa_sorted),
        groups=groups,
        has_clinical=has_clinical,
        warnings=warnings,
    )
