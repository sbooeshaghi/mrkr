"""LLM interface for mrid tool"""

import asyncio
from collections import defaultdict
from enum import Enum
from typing import Dict, List, Type, TypeVar, Union

from pandas import DataFrame
from pydantic import BaseModel, Field
from pydantic_ai import Agent

from .models import CellLabelToNameResponse, InferenceModel
from pydantic_ai.settings import ModelSettings
T = TypeVar("T", bound=BaseModel)


async def run_extraction_agent(
    model: str,
    system_prompt: str,
    user_prompt: str,
    output_type: Type[T],
) -> T:
    """
    Run an extraction agent with the given model, prompts, and output type.

    Args:
        model (str): The model to use for the agent.
        system_prompt (str): The system prompt for the agent.
        user_prompt (str): The user prompt for the agent.
        output_type (type[BaseModel]): The Pydantic model for the output.

    Returns:
        BaseModel: The Pydantic model instance containing the extracted data.
    """
    agent = Agent(model, output_type=output_type, system_prompt=system_prompt, model_settings=ModelSettings(temperature=0))
    result = await agent.run(user_prompt)
    return result.output




def extract_from_text_with_combined_enums(
    text: str,
    spec_df: "DataFrame",
    context: str = "",
    model: InferenceModel = InferenceModel.OPENAI_GPT4o,
) -> List[Dict[str, str]]:
    """Extract cell type-gene associations from text using enum-based mapping for both group_names and data_ids in a single pass."""
    return asyncio.run(
        extract_from_text_with_combined_enums_async(text, spec_df, context, model)
    )


async def extract_from_text_with_combined_enums_async(
    text: str,
    spec_df: "DataFrame",
    context: str = "",
    model: InferenceModel = InferenceModel.OPENAI_GPT4o,
) -> List[Dict[str, str]]:
    # Get unique cell type names and data_ids from spec
    unique_group_names = spec_df["group_name"].unique().tolist()
    unique_data_ids = spec_df["data_id"].unique().tolist()

    # Create tab-delimited format showing group_label -> group_name mappings
    cell_type_mappings = []
    for _, row in spec_df[["group_label", "group_name"]].drop_duplicates().iterrows():
        cell_type_mappings.append(f"{row['group_label']}\t{row['group_name']}")

    # Dynamically create enums from unique cell types and data_ids
    CellTypeSource = Enum(
        "CellTypeSource", {cell_type: cell_type for cell_type in unique_group_names}
    )

    DataIdSource = Enum("DataIdSource", {data_id: data_id for data_id in unique_data_ids})

    class DynamicMarkerGeneExtraction(BaseModel):
        group_label: str = Field(
            description="The exact cell type label as written in the text (e.g., 'T cells', 'neurons', 'hepatocytes')"
        )
        feature_label: str = Field(
            description="The exact gene name as written in the text (e.g., 'CD3D', 'GAD1', 'ALB')"
        )
        group_name: CellTypeSource = Field(  # type: ignore
            description="The enum value that best matches the cell type context from the provided cell type list"
        )
        data_id: DataIdSource = Field(  # type: ignore
            description="The enum value that best matches the data source/dataset context from the provided data_id list"
        )
        source_rationale: str = Field(
            description="The exact sentence, phrase, or contiguous text span that describes the association between the cell type and marker gene"
        )

    class DynamicMarkerGeneExtractions(BaseModel):
        extractions: List[DynamicMarkerGeneExtraction] = Field(
            description="List of extracted gene-cell type marker associations from the scientific text"
        )

    system_prompt = "You extract gene-cell type associations from scientific text. Use only the provided enum values for both cell types and data sources. Map each association to the best matching enums based on the surrounding context. Your output must reflect exactly what the text says—do not infer."

    user_prompt = f"""
You are given a scientific text. Extract all **explicit** associations between cell types and marker genes.

---

**Available Cell Types (tab-delimited: group_label → group_name):**
{chr(10).join(cell_type_mappings)}

**Available Data Sources (data_id enum):**
{chr(10).join(unique_data_ids)}

---

**Extraction Instructions**

**COMPREHENSIVE EXTRACTION**: Extract **every possible mention** of cell-type marker gene associations from the text. Be thorough and exhaustive - do not miss any gene-cell type relationships mentioned anywhere in the text.

For **each valid association**, extract:
- `feature_label`: Gene name, exactly as written in the text.
- `group_label`: Cell type name, exactly as written in the text.
- `group_name`: The best matching cell type enum value from the provided list.
- `data_id`: The best matching data source enum value from the provided list.
- `source_rationale`: The exact sentence, phrase, or text span that describes the association.

**Rules:**
1. **EXTRACT ALL MENTIONS**: Scan the entire text thoroughly. Look for gene-cell type associations in:
   - Main text paragraphs
   - Figure captions and legends
   - Table descriptions and headers
   - Methods sections
   - Results sections
   - Supplementary information references
2. Only extract associations where both a gene and a cell type are clearly mentioned in the same sentence or phrase.
3. For 1 cell type and multiple genes, or vice versa, return all valid pairs.
4. Always include the exact span of text showing the association (`source_rationale`).
5. Do not infer—only extract what is explicitly stated.

**Handling Pluralities:**
- **Multiple cell types + single gene**: If the text mentions multiple cell types expressing a single gene, create separate extraction entries for each cell type paired with that gene.
- **Single cell type + multiple genes**: If the text mentions a single cell type expressing multiple genes, create separate extraction entries for that cell type paired with each gene.
- **Multiple cell types + multiple genes**: If the text mentions multiple cell types expressing multiple genes, create extraction entries for ALL valid combinations (each cell type with each gene).
- Each extraction entry should have the same `source_rationale` (the original text span) but different `group_label`/`group_name` or `feature_label` combinations.

**Cell Type Mapping:**
The cell type list above shows mappings in tab-delimited format: group_label → group_name.
- **Important**: The text may NOT exactly match the group_label or group_name from the spec.
- Use the provided mappings as **guidance** to understand cell type relationships and select the most appropriate `group_name`.
- Always extract the `group_label` exactly as written in the text, even if it doesn't perfectly match the spec.
- For `group_name`, analyze the text context and choose the best matching value from the available group_name options (second column).
- The text might use variations, synonyms, abbreviations, or slightly different terminology than what appears in the spec.
- **Numeric group_labels**: When group_labels in the spec are simply numbers (e.g., "1", "2", "3"), interpret these as "cluster numbers" when mapping. If the text mentions "cluster 1", "cluster 2", "group 1", "subpopluation 1" etc., map accordingly to the numeric group_labels.

**General Class Label Mapping:**
Authors often refer to sets of specific cell types using broader class labels. Use these strategies to map general terms to specific cell types:
- **Biological hierarchies**: General terms like "Germ cells" may refer to multiple specific types (e.g., "SPERMATOGONIA CELL", "MEIOTIC SPERMATOCYTE", "POST-MEIOTIC ROUND SPERMATID", "ELONGATED SPERMATID"). Choose the most contextually appropriate specific type.
- **Numbered series**: Abbreviated labels like "SPG1-6" may correspond to numbered cell types in the spec (e.g., "TYPE A UNDIFFERENTIATED SPERMATOGONIAL CELL 1", "TYPE A UNDIFFERENTIATED SPERMATOGONIAL CELL 2", etc.). Match based on number and biological context.
- **Functional groups**: Terms like "immune cells", "epithelial cells", or "neural cells" should map to the most specific matching cell type from that functional category in the spec.
- **Developmental stages**: References to early/late, immature/mature, or developmental stages should map to the corresponding stage-specific cell type in the spec.
- **When multiple matches exist**: If a general term could map to several specific cell types, choose based on the surrounding context, gene associations, and biological relevance.

- Your job is to interpret the text and map it to the most semantically appropriate `group_name` from the provided options.

**Enum Mapping:**
When choosing `group_name`, match the context of the text to the cell type enum.
When choosing `data_id`, consider:
- Section titles, figure captions, and table names (e.g., "Table 2").
- Paragraphs discussing marker categories (e.g., "immune", "vascular", "global").
- The description field which provides additional context about the file or analysis.
- Broader context—e.g., if discussing immune cell types within an immune-focused section, choose appropriate data source.

**Cell Type Specificity:**
- Use the most specific label when subpopulations are mentioned.
- If the text gives a general label but enums only list subtypes, choose the best representative subtype.
- For ambiguous general terms, prefer the broader enum.

**Data Source Specificity:**
- Choose the data source that would most likely contain this specific cell type and gene association.
- Consider the context and source of the information in the text.

**Output Format:**
Your response must be a **list of JSON objects**, each matching the following schema:

```json
{{
    "group_label": "string",
    "feature_label": "string", 
    "group_name": "enum",
    "data_id": "enum",
    "source_rationale": "string"
}}
```

The list must be wrapped inside a single top-level field called `extractions`. For example:

```json
{{
    "extractions": [
        {{
            "group_label": "T cells",
            "feature_label": "CD3D",
            "group_name": "T cell",
            "data_id": "immune_data",
            "source_rationale": "CD3D expression was found to be enriched in T cells."
        }}
    ]
}}
```

**Important:**
- Stick to what's written. No guessing.
- Your job is extraction and precise mapping to both cell type and data source enums.
- If unsure, pick the closest enum based on context.
- Extract all possible cell-type marker gene pairs!

{f"Additional Context: {context}" if context else ""}

Text to Analyze:
{text}
"""

    try:
        # Note: The model parameter is now a string, e.g., "openai:gpt-4o"
        extraction_result = await run_extraction_agent(
            model=model.value,
            system_prompt=system_prompt,
            user_prompt=user_prompt,
            output_type=DynamicMarkerGeneExtractions,
        )
        # The agent returns a Pydantic model, so we extract the list of dicts.
        return [item.model_dump() for item in extraction_result.extractions]
    except Exception as e:
        print(f"Warning: Error extracting from text with combined enums: {e}")
        return []


def infer_cell_names(
    cell_labels: List[str], model: InferenceModel, context: str = ""
) -> List[Dict[str, Union[str, None]]]:
    """
    Infer cell names for a list of cell labels using language models.

    Args:
        cell_labels: List of cell labels to infer from
        model: The model to use for inference
        context: Additional context to help with inference

    Returns:
        List of dictionaries with original_label and standardized_name
    """
    # Convert to set for unique labels
    unique_labels = set(cell_labels)

    # Filter out generic cluster names
    labels_to_infer = {
        label
        for label in unique_labels
        if not any(
            label.lower().startswith(prefix)
            for prefix in ["cluster", "group", "cell", "type"]
        )
    }
    print(labels_to_infer)

    if not labels_to_infer:
        return [
            {"original_label": label, "standardized_name": label, "ontology_id": None}
            for label in cell_labels
        ]

    system_prompt = """You are an expert at standardizing cell type names in single-cell genomics.
Given a list of cell labels, you will return standardized cell names.
If you cannot determine a meaningful cell type name, return the original label.
If the cell type label is an acronym, attempt to infer the name from the acronym.
If the cell type label as a number associated with it, include it in the cell type name.
Always maintain the exact same number of entries in your response as in the input.
You must respond with a valid JSON object containing a list of mappings between original labels and inferred names.
The same input label may correspond to multiple inferred labels. Return all of them as separate entries.
"""

    user_prompt = f"""Given the following cell labels, return standardized cell names for each.
If you cannot determine a meaningful cell type name, return the original label.

Examples:
{{
        "mappings": [
        {{"original_label": "CD4+ T cells", "inferred_name": "CD4 T cell"}},
        {{"original_label": "B lymphocytes", "inferred_name": "B cell"}},
        {{"original_label": "cluster1", "inferred_name": "cluster1"}},
        {{
            "original_label": "Amacrine/Horizontal cells", "inferred_name": "Amacrine cells"}},
        {{
            "original_label": "Amacrine/Horizontal cells", "inferred_name": "Horizontal cells"}}
    ]
}}

{f"Additional context: {context}" if context else ""}

Input labels:
{chr(10).join(labels_to_infer)}

Return a JSON object with this structure:
{{
    "mappings": [
        {{"original_label": "label1", "inferred_name": "name1"}},
        {{"original_label": "label2/label3", "inferred_name": "name2"}}
        {{"original_label": "label2/label3", "inferred_name": "name3"}}
    ]
}}"""

    try:
        extraction_result = asyncio.run(
            run_extraction_agent(
                model=model.value,
                system_prompt=system_prompt,
                user_prompt=user_prompt,
                output_type=CellLabelToNameResponse,
            )
        )
        # Create mapping for inferred labels
        inferred_names: Dict[str, List[str]] = defaultdict(list)
        for mapping in extraction_result.mappings:
            if (
                any(
                    mapping.inferred_name.lower().startswith(prefix)
                    for prefix in ["cluster", "group", "cell", "type"]
                )
                or mapping.inferred_name.lower()
                == mapping.original_label.lower()
            ):
                inferred_names[mapping.original_label].append(
                    mapping.original_label
                )
            else:
                inferred_names[mapping.original_label].append(
                    mapping.inferred_name
                )

        # Add original labels for generic cluster names
        for label in unique_labels:
            if label not in inferred_names:
                inferred_names[label].append(label)

        # Map back to original list order and take first inferred name
        result = []
        for label in cell_labels:
            if label in inferred_names and inferred_names[label]:
                standardized_name = inferred_names[label][
                    0
                ]  # Take first inferred name
            else:
                standardized_name = label

            result.append(
                {
                    "original_label": label,
                    "standardized_name": standardized_name,
                    "ontology_id": None,  # Keep for backward compatibility
                }
            )

        return result

    except Exception as e:
        print(f"Warning: Error during inference: {e}")
        return [
            {"original_label": label, "standardized_name": label, "ontology_id": None}
            for label in cell_labels
        ]
