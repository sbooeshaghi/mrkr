"""Typed language-model output for marker claims."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field, model_validator


class TermExtraction(BaseModel):
    """One biological term emitted by the extraction model."""

    sub_span: str | None = Field(
        default=None,
        description="Exact term text within span_literal, or null when implicit.",
    )
    normalized_label: str
    term_type: Literal["gene", "celltype", "comparison", "tissue"]
    direction: Literal["positive", "negative"] | None = None

    @model_validator(mode="after")
    def validate_direction(self):
        if self.term_type == "gene" and self.direction is None:
            raise ValueError("gene terms require positive or negative direction")
        if self.term_type != "gene" and self.direction is not None:
            raise ValueError("direction is only valid for gene terms")
        return self


class ClaimExtraction(BaseModel):
    """One target cell type and its marker evidence."""

    span_literal: str
    summary: str
    terms: list[TermExtraction]


class ClaimsResult(BaseModel):
    """Complete model response for one manuscript."""

    claims: list[ClaimExtraction]
