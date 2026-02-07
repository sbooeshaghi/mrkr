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
    source_type: Literal["text", "image", "deg"]
    source_rationale: str
    source_id: str
    data_id: Optional[str] = None

    # DEG-specific metrics (null for text/image sources)
    metrics_pcorr: Optional[float] = None
    metrics_logfc: Optional[float] = None
    metrics_rank: Optional[int] = None
