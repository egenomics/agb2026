"""Pydantic data models shared by parsers.py and the charts package."""
from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field

# ── Raw upload payloads ───────────────────────────────────────────────────────

class TSVUpload(BaseModel):
    """A raw TSV file sent from the browser as a string."""

    content: str = Field(..., description="Full TSV file content as UTF-8 string")


# ── Parsed data shapes ────────────────────────────────────────────────────────

class FeatureTableResult(BaseModel):
    """Output of parsing a QIIME2 feature-table.tsv.

    counts[feature_id][sample_id] = raw count.
    """

    samples:    list[str]                   = Field(..., description="Ordered sample IDs")
    features:   list[str]                   = Field(..., description="Ordered ASV/OTU IDs")
    counts:     dict[str, dict[str, float]] = Field(..., description="counts[feature_id][sample_id]")
    n_samples:  int                         = Field(..., description="Number of samples")
    n_features: int                         = Field(..., description="Number of features (ASVs/OTUs)")


class TaxonomyResult(BaseModel):
    """Output of parsing a QIIME2 taxonomy.tsv.

    assignments[feature_id] = family-level name.
    """

    assignments:      dict[str, str] = Field(..., description="feature_id -> family name")
    unclassified_pct: float          = Field(..., description="% features unclassified at family level")
    unique_families:  list[str]      = Field(..., description="Sorted unique family names found")


class SampleMetadata(BaseModel):
    """Metadata for a single sample."""

    sample_id:  str
    group:      str        = Field(..., description="Full subgroup label e.g. EAA_T0")
    base_group: str        = Field(..., description="Treatment group e.g. EAA")
    timepoint:  str        = Field(..., description="Timepoint label e.g. T0, T84")
    time:       int | None = Field(None, description="Numeric time in days, None if not parseable")
    patient:    str        = Field(..., description="Patient/subject identifier")
    sixmwt:     float      = Field(0.0, description="6-min walk test (m), 0 if absent")
    il18:       float      = Field(0.0, description="IL-18 cytokine (pg/mL), 0 if absent")


class MetadataResult(BaseModel):
    """Output of parsing a QIIME2 metadata.tsv."""

    samples:      list[SampleMetadata]
    groups:       list[str] = Field(..., description="Unique group values found")
    base_groups:  list[str] = Field(..., description="Unique base_group values found")
    timepoints:   list[str] = Field(..., description="Unique timepoint values found")
    has_clinical: bool      = Field(..., description="True if sixmwt or il18 data found")
    n_samples:    int


class AlphaDiversityEntry(BaseModel):
    """Alpha diversity metrics for a single sample."""

    sample_id: str
    shannon:   float = 0.0
    simpson:   float = 0.0
    observed:  float = 0.0
    faith_pd:  float = 0.0
    pielou:    float = 0.0


class AlphaDiversityResult(BaseModel):
    """Output of parsing a QIIME2 alpha-diversity.tsv."""

    samples:         list[AlphaDiversityEntry]
    metrics_present: list[str] = Field(..., description="Columns found in the file")


class DistanceMatrixResult(BaseModel):
    """Output of parsing a QIIME2 distance-matrix.tsv."""

    samples: list[str]
    matrix:  list[list[float]] = Field(..., description="Symmetric n x n distance matrix")
    n:       int                = Field(..., description="Number of samples")


# ── Integrated sample row (what charts consume) ───────────────────────────────

class SampleRow(BaseModel):
    """A fully integrated, chart-ready row.

    Produced by parsers.integrate() — combines feature table, taxonomy,
    metadata, and optional alpha diversity.  Taxon abundances are stored
    as dynamic extra fields.
    """

    model_config = ConfigDict(extra="allow")

    sample_id:  str
    patient:    str
    group:      str
    base_group: str
    timepoint:  str
    time:       int | None
    shannon:    float = 0.0
    simpson:    float = 0.0
    pielou:     float = 0.0
    observed:   float = 0.0
    faith_pd:   float = 0.0
    sixmwt:     float = 0.0
    il18:       float = 0.0


class IntegrateResult(BaseModel):
    """Output of parsers.integrate() — ready for all chart builders."""

    rows:         list[SampleRow]
    taxa:         list[str]  = Field(..., description="Family names sorted by mean abundance desc")
    n_samples:    int
    n_taxa:       int
    groups:       list[str]
    has_clinical: bool
    warnings:     list[str] = Field(default_factory=list)


# ── Validation errors ─────────────────────────────────────────────────────────

class ParseError(BaseModel):
    """Structured error returned when a parse fails."""

    file_type: str
    message:   str
    line:      int | None = None
    hint:      str | None = None
