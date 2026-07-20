"""Data models for mrkr."""

from typing import Optional, Literal, List
from pydantic import BaseModel, Field


class ExtractionResponse(BaseModel):
    """Single extraction from LLM (minimal fields)."""

    organism: str = Field(
        description="Species Latin name (e.g., homo_sapiens, mus_musculus)"
    )
    group_label: str = Field(
        description="Exact cell type text as it appears in source"
    )
    group_name: str = Field(
        description="Normalized cell type name (UPPERCASE)"
    )
    feature_label: str = Field(
        description="Exact gene name as it appears in source"
    )
    feature_name: str = Field(
        description="Normalized gene name (UPPERCASE)"
    )
    source_type: Literal["text", "image"] = Field(
        description="Where this was extracted from (text or image)"
    )
    source_rationale: str = Field(
        description="Text snippet or figure description showing this relationship"
    )
    data_id: Optional[str] = Field(
        default=None,
        description="Data source identifier (e.g., file#sheet from DEG tables)"
    )


class ExtractionsResult(BaseModel):
    """Container for all extractions from LLM."""

    extractions: List[ExtractionResponse] = Field(
        description="List of all extracted marker gene associations"
    )


# --- Claim extraction models (the current format; LLM output = spans + labels, NO ids) ---


class TermExtraction(BaseModel):
    """A single grounded-claim term as emitted by the LLM. No ontology id — grounding adds it."""

    sub_span: Optional[str] = Field(
        default=None,
        description="Verbatim substring of this claim's span_literal (the surface token), "
                    "or null if the entity is not present in the sentence (implicit)",
    )
    normalized_label: str = Field(
        description="Canonical name; MUST be a verbatim substring of this claim's summary"
    )
    term_type: Literal["gene", "celltype", "comparison", "tissue"] = Field(
        description="gene | celltype | comparison | tissue"
    )
    direction: Optional[Literal["positive", "negative"]] = Field(
        default=None,
        description="Gene terms only: 'negative' if the cell type does NOT express the gene, "
                    "else 'positive'",
    )


class ClaimExtraction(BaseModel):
    """A single marker claim as emitted by the LLM (one target cell type + its marker terms)."""

    span_literal: str = Field(
        description="Verbatim, exact, contiguous substring of the paper stating the marker(s)"
    )
    summary: str = Field(
        description="Normalized rewrite of span_literal; each normalized_label appears in it verbatim"
    )
    terms: List[TermExtraction] = Field(
        description="Exactly one celltype term (the target) + one or more gene terms + "
                    "optional comparison/tissue terms"
    )


class ClaimsResult(BaseModel):
    """Container for all marker claims from the LLM."""

    claims: List[ClaimExtraction] = Field(
        description="List of marker claim objects"
    )


class Evidence(BaseModel):
    """Final output format for evidence record (uniform structure for all sources)."""

    organism: str
    group_label: str
    group_name: str
    group_id: Optional[str] = None
    feature_label: str
    feature_name: str
    feature_id: Optional[str] = None
    source_type: Literal["text", "image", "deg", "generated", "predicted", "selected"]
    source_rationale: str
    source_id: str
    data_id: Optional[str] = None

    # DEG-specific metrics (null for text/image sources)
    metrics_pcorr: Optional[float] = None
    metrics_logfc: Optional[float] = None
    metrics_rank: Optional[int] = None


# --- Generation models (for mrkr generate) ---


class GeneratedMarker(BaseModel):
    """A single generated marker gene association."""

    group_name: str = Field(description="Cell type name (UPPERCASE)")
    feature_name: str = Field(description="Gene symbol (UPPERCASE)")
    rationale: str = Field(
        description="Brief explanation of why this gene is a marker for this cell type"
    )


class GenerationResult(BaseModel):
    """Container for generated marker genes (celltypes-to-genes mode)."""

    generations: List[GeneratedMarker] = Field(
        description="List of generated marker gene associations"
    )


class CellTypePrediction(BaseModel):
    """A single cell type prediction for a gene group."""

    group_id: str = Field(description="The anonymous group identifier (e.g., 'Group 1')")
    predicted_cell_type: str = Field(description="Predicted cell type name (UPPERCASE)")
    rationale: str = Field(
        description="Brief explanation of why these genes indicate this cell type"
    )


class CellTypePredictionResult(BaseModel):
    """Container for cell type predictions (genes-to-celltypes mode)."""

    predictions: List[CellTypePrediction] = Field(
        description="List of cell type predictions"
    )


# --- Selection models (for mrkr select) ---


class SelectedMarker(BaseModel):
    """A marker gene selected from DEG data."""

    group_name: str = Field(
        description="Cell type name (UPPERCASE, exactly as in input)"
    )
    feature_name: str = Field(
        description="Gene symbol (UPPERCASE, exactly as in input)"
    )
    rationale: str = Field(
        description="Brief explanation of why this gene is a marker for this cell type"
    )


class SelectionResult(BaseModel):
    """Container for selected marker genes."""

    selections: List[SelectedMarker] = Field(
        description="List of selected markers"
    )
