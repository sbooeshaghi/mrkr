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


class Evidence(BaseModel):
    """Final output format for evidence record (uniform structure for all sources)."""

    organism: str
    group_label: str
    group_name: str
    group_id: Optional[str] = None
    feature_label: str
    feature_name: str
    feature_id: Optional[str] = None
    source_type: Literal["text", "image", "deg", "generated", "predicted"]
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
