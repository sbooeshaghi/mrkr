"""Generation functionality for mrkr generate command."""

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Type

import anthropic
from anthropic.types import Message
from pydantic import BaseModel

from .config import config
from .llm import load_prompt_template, extract_json_from_response
from .metrics import save_metrics, Timer
from .models import (
    GenerationResult,
    CellTypePredictionResult,
)


def load_input_cell_types(path: Path) -> List[str]:
    """Extract unique group_names from an extracted.json file.

    Args:
        path: Path to extracted.json (list of evidence records)

    Returns:
        Sorted list of unique group_name values
    """
    data = json.loads(path.read_text(encoding="utf-8"))
    group_names = set()
    for record in data:
        gn = record.get("group_name")
        if gn:
            gn = gn.strip().upper()
            if gn:
                group_names.add(gn)
    return sorted(group_names)


def load_input_gene_groups(path: Path) -> Dict[str, List[str]]:
    """Build gene groups from an extracted.json file.

    Args:
        path: Path to extracted.json (list of evidence records)

    Returns:
        Dict mapping group_name → sorted list of unique feature_names
    """
    data = json.loads(path.read_text(encoding="utf-8"))
    groups: Dict[str, set] = {}
    for record in data:
        gn = (record.get("group_name") or "").strip().upper()
        fn = (record.get("feature_name") or "").strip().upper()
        if gn and fn:
            if gn not in groups:
                groups[gn] = set()
            groups[gn].add(fn)
    return {gn: sorted(fns) for gn, fns in sorted(groups.items())}


def call_claude_generate(
    prompt: str,
    response_model: Type[BaseModel],
    verbose: bool = False,
) -> Tuple[BaseModel, Message]:
    """Call Claude API expecting a JSON response matching the given Pydantic model.

    Args:
        prompt: The prompt text
        response_model: Pydantic model class for response validation
        verbose: Whether to print progress

    Returns:
        Tuple of (parsed response model instance, Message)
    """
    client = anthropic.Anthropic(api_key=config.anthropic_api_key)

    if verbose:
        print(f"   Calling Claude ({config.anthropic_model})...")

    with client.messages.stream(
        model=config.anthropic_model,
        max_tokens=32000,
        temperature=0.0,
        messages=[{"role": "user", "content": prompt}],
    ) as stream:
        message = stream.get_final_message()

    response_text = message.content[0].text

    if verbose:
        print(
            f"   Response: {message.usage.input_tokens} tokens in, "
            f"{message.usage.output_tokens} tokens out"
        )

    # Parse JSON response
    json_text = extract_json_from_response(response_text)

    if not json_text or not json_text.strip():
        raise ValueError(
            f"Failed to extract JSON from response. Response text: {response_text[:1000]}"
        )

    try:
        result = response_model.model_validate_json(json_text)
        return result, message
    except Exception:
        pass

    # Fallback: try json_repair
    from json_repair import repair_json

    repaired = repair_json(json_text, return_objects=True)
    if verbose:
        print(f"   JSON was malformed, repaired via json_repair")
    result = response_model.model_validate(repaired)
    return result, message


def generate_celltypes_to_genes(
    input_path: Path,
    species: str,
    verbose: bool = False,
    metrics_path: Optional[Path] = None,
    command: str = "",
) -> List[dict]:
    """Generate marker genes for cell types using LLM training knowledge.

    Args:
        input_path: Path to extracted.json with group_name fields
        species: Species Latin name (e.g., homo_sapiens)
        verbose: Whether to print progress
        metrics_path: Optional path to save LLM metrics
        command: CLI command string for metrics

    Returns:
        List of evidence dictionaries with source_type="generated"
    """
    cell_types = load_input_cell_types(input_path)

    if verbose:
        print(f"\n   Loaded {len(cell_types)} cell types from {input_path.name}")
        for ct in cell_types:
            print(f"     - {ct}")

    # Format cell types for prompt
    cell_types_text = "\n".join(f"- {ct}" for ct in cell_types)

    # Build prompt
    template = load_prompt_template("generate_celltypes_to_genes")
    prompt = template.format(species=species, cell_types=cell_types_text)

    # Call Claude
    with Timer() as timer:
        result, message = call_claude_generate(
            prompt, GenerationResult, verbose=verbose
        )

    # Save metrics if requested
    if metrics_path:
        save_metrics(
            output_path=metrics_path,
            model=config.anthropic_model,
            message=message,
            processing_time_sec=timer.elapsed,
            num_extractions=len(result.generations),
            command=command,
        )
        if verbose:
            print(f"   Saved metrics to: {metrics_path}")

    # Convert to evidence dicts
    source_id = f"generated:{config.anthropic_model}"
    records = []
    for gen in result.generations:
        records.append({
            "organism": species,
            "group_label": gen.group_name.strip().upper(),
            "group_name": gen.group_name.strip().upper(),
            "group_id": None,
            "feature_label": gen.feature_name.strip().upper(),
            "feature_name": gen.feature_name.strip().upper(),
            "feature_id": None,
            "source_type": "generated",
            "source_rationale": gen.rationale.strip(),
            "source_id": source_id,
            "data_id": None,
            "metrics_pcorr": None,
            "metrics_logfc": None,
            "metrics_rank": None,
        })

    if verbose:
        print(f"   Generated {len(records)} marker gene associations")

    return records


def generate_genes_to_celltypes(
    input_path: Path,
    species: str,
    verbose: bool = False,
    metrics_path: Optional[Path] = None,
    command: str = "",
) -> List[dict]:
    """Predict cell types from gene groups using LLM training knowledge.

    Groups are anonymized (Group 1, Group 2, ...) so the LLM must predict
    from gene content alone.

    Args:
        input_path: Path to extracted.json with group_name/feature_name fields
        species: Species Latin name (e.g., homo_sapiens)
        verbose: Whether to print progress
        metrics_path: Optional path to save LLM metrics
        command: CLI command string for metrics

    Returns:
        List of evidence dictionaries with source_type="predicted"
    """
    gene_groups = load_input_gene_groups(input_path)

    if verbose:
        print(f"\n   Loaded {len(gene_groups)} gene groups from {input_path.name}")

    # Anonymize: map real group names to "Group N"
    group_names_ordered = list(gene_groups.keys())
    anon_map = {}  # "Group N" → real group_name
    anon_lines = []
    for i, gn in enumerate(group_names_ordered, 1):
        anon_id = f"Group {i}"
        anon_map[anon_id] = gn
        genes = ", ".join(gene_groups[gn])
        anon_lines.append(f"- {anon_id}: {genes}")
        if verbose:
            print(f"     {anon_id} ({gn}): {len(gene_groups[gn])} genes")

    gene_groups_text = "\n".join(anon_lines)

    # Build prompt
    template = load_prompt_template("generate_genes_to_celltypes")
    prompt = template.format(species=species, gene_groups=gene_groups_text)

    # Call Claude
    with Timer() as timer:
        result, message = call_claude_generate(
            prompt, CellTypePredictionResult, verbose=verbose
        )

    # Save metrics if requested
    if metrics_path:
        save_metrics(
            output_path=metrics_path,
            model=config.anthropic_model,
            message=message,
            processing_time_sec=timer.elapsed,
            num_extractions=len(result.predictions),
            command=command,
        )
        if verbose:
            print(f"   Saved metrics to: {metrics_path}")

    # Convert to evidence dicts
    source_id = f"predicted:{config.anthropic_model}"
    records = []
    for pred in result.predictions:
        group_id = pred.group_id.strip()
        real_group_name = anon_map.get(group_id, group_id)
        genes = gene_groups.get(real_group_name, [])

        # One record per gene in the group, with the predicted cell type
        for gene in genes:
            records.append({
                "organism": species,
                "group_label": pred.predicted_cell_type.strip().upper(),
                "group_name": pred.predicted_cell_type.strip().upper(),
                "group_id": None,
                "feature_label": gene,
                "feature_name": gene,
                "feature_id": None,
                "source_type": "predicted",
                "source_rationale": pred.rationale.strip(),
                "source_id": source_id,
                "data_id": None,
                "metrics_pcorr": None,
                "metrics_logfc": None,
                "metrics_rank": None,
                "_original_group_name": real_group_name,
            })

    if verbose:
        # Summarize predictions
        for pred in result.predictions:
            group_id = pred.group_id.strip()
            real_name = anon_map.get(group_id, "?")
            print(f"     {group_id} ({real_name}) -> {pred.predicted_cell_type}")
        print(f"   Generated {len(records)} predicted marker gene associations")

    return records
