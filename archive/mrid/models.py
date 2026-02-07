"""Data models for mrid tool"""

from enum import Enum
from typing import Dict, List, Optional

from pydantic import BaseModel, Field


class InferenceModel(Enum):
    """Available LLM models for inference"""

    OPENAI_GPT3_5_TURBO = "openai:gpt-3.5-turbo"
    OPENAI_GPT4 = "openai:gpt-4"
    OPENAI_GPT4o = "openai:gpt-4o"
    OPENAI_GPT4_1 = "openai:gpt-4.1"
    OLLAMA_LLAMA3 = "ollama:llama3"
    GEMINI_1_5_FLASH = "google-gla:gemini-1.5-flash"


class CellLabelToName(BaseModel):
    """Structured output for cell name inference."""

    original_label: str = Field(
        ..., description="The original cell label from the data"
    )
    inferred_name: str = Field(
        ..., description="The standardized cell name inferred from the label"
    )


class CellLabelToNameResponse(BaseModel):
    """Response format for cell name inference."""

    mappings: List[CellLabelToName] = Field(
        ..., description="List of cell label to name mappings"
    )


class CellTypeGeneAssociation(BaseModel):
    """Structured output for cell type-gene associations"""

    cell_type: str
    gene: str
    context: str
    data_id: str


class CellTypeGeneAssociationResponse(BaseModel):
    """Response containing multiple cell type-gene associations"""

    associations: List[CellTypeGeneAssociation]


class GeneMapping(BaseModel):
    """Structured output for gene name mapping"""

    original_name: str
    standardized_name: Optional[str] = None
    gene_id: Optional[str] = None


class GeneMappingResponse(BaseModel):
    """Response containing multiple gene mappings"""

    genes: List[GeneMapping]


class MarkerGeneExtraction(BaseModel):
    """Structured output for marker gene extraction with enum-based cell type mapping"""

    group_label: str = Field(
        description="The exact cell type label as written in the text (e.g., 'T cells', 'neurons', 'hepatocytes')"
    )
    feature_label: str = Field(
        description="The exact gene name as written in the text (e.g., 'CD3D', 'GAD1', 'ALB')"
    )
    group_name: str = Field(
        description="The best matching cell type name from the provided list"
    )
    source_rationale: str = Field(
        description="The exact sentence, phrase, or contiguous text span that describes the association between the cell type and marker gene"
    )


class MarkerGeneExtractions(BaseModel):
    """Response containing multiple marker gene extractions"""

    extractions: List[MarkerGeneExtraction] = Field(
        description="List of extracted gene-cell type marker associations from the scientific text"
    )


class Evidence(BaseModel):
    """Evidence record for gene-cell type associations"""

    organism: str
    group_label: str
    feature_label: str
    group_name: Optional[str] = None
    group_id: Optional[str] = None
    feature_name: Optional[str] = None
    feature_id: Optional[str] = None
    source_type: str = "deg"
    source_rationale: str = "extracted from deg"
    source_id: str = ""
    data_id: str = ""
    metrics_pcorr: Optional[float] = None
    metrics_logfc: Optional[float] = None
    metrics_rank: Optional[int] = None


class SpecRecord(BaseModel):
    """Specification record for file metadata"""

    file_id: str
    file_name: str
    file_type: str
    file_uri: str
    data_id: str
    data_type: str
    group_label: str
    group_name: str
    group_id: str
