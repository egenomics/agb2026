"""
app/routes/parse.py

All /api/parse/* endpoints.
Each endpoint:
  - Accepts a TSV file as multipart/form-data OR as raw string in JSON body
  - Validates input
  - Delegates ALL logic to services/parsers.py
  - Returns typed response or structured error

No business logic lives here.
"""

from __future__ import annotations

import logging
from typing import Optional

from fastapi import APIRouter, File, HTTPException, UploadFile, status

from app.models.sample import (
    FeatureTableResult,
    TaxonomyResult,
    MetadataResult,
    AlphaDiversityResult,
    DistanceMatrixResult,
    IntegrateResult,
    ParseError,
)
from app.services import parsers

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/api/parse", tags=["parse"])


# ── Helper: read upload as UTF-8 string ──────────────────────────────────────

async def _read_upload(file: UploadFile) -> str:
    raw = await file.read()
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


def _parse_error(file_type: str, exc: ValueError) -> HTTPException:
    """Convert a service-layer ValueError into a 422 HTTPException with context."""
    return HTTPException(
        status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
        detail=ParseError(
            file_type=file_type,
            message=str(exc),
            hint=_hints.get(file_type),
        ).model_dump(),
    )


_hints: dict[str, str] = {
    "feature-table": (
        "Export from QIIME2 with: "
        "qiime tools export --input-path table.qza --output-path . "
        "followed by biom convert -i feature-table.biom -o feature-table.tsv --to-tsv"
    ),
    "taxonomy": (
        "Export from QIIME2 with: "
        "qiime tools export --input-path taxonomy.qza --output-path ."
    ),
    "metadata": (
        "Metadata must have a 'sample-id' column. "
        "See: https://docs.qiime2.org/2024.5/tutorials/metadata/"
    ),
    "alpha-diversity": (
        "Export from QIIME2 with: "
        "qiime tools export --input-path shannon_vector.qza --output-path ."
    ),
    "distance-matrix": (
        "Export from QIIME2 with: "
        "qiime tools export --input-path bray_curtis_distance_matrix.qza --output-path ."
    ),
}


# ── Feature table ─────────────────────────────────────────────────────────────

@router.post(
    "/feature-table",
    response_model=FeatureTableResult,
    summary="Parse a QIIME2 feature-table.tsv",
    description=(
        "Accepts a QIIME2-exported feature-table.tsv (converted from BIOM). "
        "Returns sample IDs, feature IDs, and a count matrix."
    ),
)
async def parse_feature_table(
    file: UploadFile = File(..., description="feature-table.tsv file"),
) -> FeatureTableResult:
    content = await _read_upload(file)
    try:
        return parsers.parse_feature_table(content)
    except ValueError as exc:
        raise _parse_error("feature-table", exc)


# ── Taxonomy ──────────────────────────────────────────────────────────────────

@router.post(
    "/taxonomy",
    response_model=TaxonomyResult,
    summary="Parse a QIIME2 taxonomy.tsv",
    description=(
        "Accepts a QIIME2 taxonomy.tsv. "
        "Extracts family-level classification for each feature. "
        "Returns feature_id → family mapping plus summary statistics."
    ),
)
async def parse_taxonomy(
    file: UploadFile = File(..., description="taxonomy.tsv file"),
) -> TaxonomyResult:
    content = await _read_upload(file)
    try:
        return parsers.parse_taxonomy(content)
    except ValueError as exc:
        raise _parse_error("taxonomy", exc)


# ── Metadata ──────────────────────────────────────────────────────────────────

@router.post(
    "/metadata",
    response_model=MetadataResult,
    summary="Parse a QIIME2 sample metadata.tsv",
    description=(
        "Accepts a QIIME2 sample metadata file. "
        "Detects group, timepoint, patient, and clinical columns automatically. "
        "Returns typed SampleMetadata for each sample."
    ),
)
async def parse_metadata(
    file: UploadFile = File(..., description="metadata.tsv file"),
) -> MetadataResult:
    content = await _read_upload(file)
    try:
        return parsers.parse_metadata(content)
    except ValueError as exc:
        raise _parse_error("metadata", exc)


# ── Alpha diversity ───────────────────────────────────────────────────────────

@router.post(
    "/alpha-diversity",
    response_model=AlphaDiversityResult,
    summary="Parse a QIIME2 alpha-diversity.tsv",
    description=(
        "Accepts a QIIME2 alpha diversity export. "
        "Detects shannon, simpson, observed_features, faith_pd, pielou columns. "
        "Multiple metrics can be present in a single file."
    ),
)
async def parse_alpha_diversity(
    file: UploadFile = File(..., description="alpha-diversity.tsv file"),
) -> AlphaDiversityResult:
    content = await _read_upload(file)
    try:
        return parsers.parse_alpha_diversity(content)
    except ValueError as exc:
        raise _parse_error("alpha-diversity", exc)


# ── Distance matrix ───────────────────────────────────────────────────────────

@router.post(
    "/distance-matrix",
    response_model=DistanceMatrixResult,
    summary="Parse a QIIME2 distance-matrix.tsv",
    description=(
        "Accepts a QIIME2 distance matrix export (e.g. Bray-Curtis, UniFrac). "
        "Returns sample IDs and the n×n distance matrix."
    ),
)
async def parse_distance_matrix(
    file: UploadFile = File(..., description="distance-matrix.tsv file"),
) -> DistanceMatrixResult:
    content = await _read_upload(file)
    try:
        return parsers.parse_distance_matrix(content)
    except ValueError as exc:
        raise _parse_error("distance-matrix", exc)


# ── Integrate (the main endpoint MicroSee.html calls) ────────────────────────

@router.post(
    "/integrate",
    response_model=IntegrateResult,
    summary="Integrate all QIIME2 files into chart-ready data",
    description=(
        "The primary endpoint. "
        "Upload feature-table, taxonomy, and metadata (required). "
        "Alpha-diversity and distance-matrix are optional. "
        "Returns a list of SampleRows with taxon abundances and diversity metrics, "
        "ready to be passed directly to chart endpoints."
    ),
)
async def integrate(
    feature_table:   UploadFile           = File(...,  description="feature-table.tsv"),
    taxonomy:        UploadFile           = File(...,  description="taxonomy.tsv"),
    metadata:        UploadFile           = File(...,  description="metadata.tsv"),
    alpha_diversity: Optional[UploadFile] = File(None, description="alpha-diversity.tsv (optional)"),
    distance_matrix: Optional[UploadFile] = File(None, description="distance-matrix.tsv (optional)"),
) -> IntegrateResult:

    # Parse each file, surfacing per-file errors clearly
    try:
        feat = parsers.parse_feature_table(await _read_upload(feature_table))
    except ValueError as exc:
        raise _parse_error("feature-table", exc)

    try:
        tax = parsers.parse_taxonomy(await _read_upload(taxonomy))
    except ValueError as exc:
        raise _parse_error("taxonomy", exc)

    try:
        meta = parsers.parse_metadata(await _read_upload(metadata))
    except ValueError as exc:
        raise _parse_error("metadata", exc)

    alpha = None
    if alpha_diversity:
        try:
            alpha = parsers.parse_alpha_diversity(await _read_upload(alpha_diversity))
        except ValueError as exc:
            raise _parse_error("alpha-diversity", exc)

    # Integrate
    try:
        return parsers.integrate(feat, tax, meta, alpha)
    except ValueError as exc:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=ParseError(
                file_type="integration",
                message=str(exc),
                hint="Ensure sample IDs match exactly between feature-table.tsv and metadata.tsv",
            ).model_dump(),
        )
