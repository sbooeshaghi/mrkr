"""Selection functionality for mrkr select command."""

import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .config import config
from .generate import call_claude_generate
from .llm import load_prompt_template
from .metrics import save_metrics, Timer
from .models import SelectionResult


def prepare_deg_data_for_prompt(
    deg_records: List[dict],
    top_n: int,
) -> Tuple[str, Dict[str, List[dict]]]:
    """Format DEG data for the selection prompt.

    Groups records by group_name, takes top N by rank (lowest rank = highest LFC),
    and formats as tab-separated sections.

    Args:
        deg_records: List of evidence dicts from DEG extraction
        top_n: Number of top DEGs per cell type to include

    Returns:
        Tuple of (formatted text for prompt, dict of group_name -> filtered records)
    """
    # Group by group_name
    by_group: Dict[str, List[dict]] = defaultdict(list)
    for rec in deg_records:
        gn = (rec.get("group_name") or "").strip().upper()
        if gn and rec.get("metrics_rank") is not None:
            by_group[gn].append(rec)

    # Sort each group by rank (ascending) and take top N
    filtered: Dict[str, List[dict]] = {}
    sections = []

    for gn in sorted(by_group.keys()):
        records = sorted(by_group[gn], key=lambda r: r["metrics_rank"])
        top_records = records[:top_n]
        filtered[gn] = top_records

        # Format as tab-separated table
        lines = [f"## {gn}"]
        lines.append("gene\tlogFC\tp_adj")
        for rec in top_records:
            gene = rec.get("feature_name", "")
            logfc = rec.get("metrics_logfc")
            pcorr = rec.get("metrics_pcorr")
            logfc_str = f"{logfc:.4f}" if logfc is not None else "NA"
            pcorr_str = f"{pcorr:.2e}" if pcorr is not None else "NA"
            lines.append(f"{gene}\t{logfc_str}\t{pcorr_str}")

        sections.append("\n".join(lines))

    return "\n\n".join(sections), filtered


def _format_deg_sections(filtered_records: Dict[str, List[dict]], cell_types: List[str]) -> str:
    """Format DEG data for a subset of cell types."""
    sections = []
    for gn in cell_types:
        recs = filtered_records[gn]
        lines = [f"## {gn}"]
        lines.append("gene\tlogFC\tp_adj")
        for rec in recs:
            gene = rec.get("feature_name", "")
            logfc = rec.get("metrics_logfc")
            pcorr = rec.get("metrics_pcorr")
            logfc_str = f"{logfc:.4f}" if logfc is not None else "NA"
            pcorr_str = f"{pcorr:.2e}" if pcorr is not None else "NA"
            lines.append(f"{gene}\t{logfc_str}\t{pcorr_str}")
        sections.append("\n".join(lines))
    return "\n\n".join(sections)


def _estimate_tokens(filtered_records: Dict[str, List[dict]], cell_types: List[str]) -> int:
    """Estimate token count for a set of cell types (~18 tokens per row + overhead)."""
    total_rows = sum(len(filtered_records[gn]) for gn in cell_types)
    return total_rows * 18 + 1000


def _build_batches(filtered_records: Dict[str, List[dict]], max_tokens: int = 160000) -> List[List[str]]:
    """Split cell types into batches that fit within token budget."""
    all_cts = sorted(filtered_records.keys())
    batches = []
    current_batch = []
    current_tokens = 1000  # prompt overhead

    for gn in all_cts:
        ct_tokens = len(filtered_records[gn]) * 18
        if current_batch and current_tokens + ct_tokens > max_tokens:
            batches.append(current_batch)
            current_batch = [gn]
            current_tokens = 1000 + ct_tokens
        else:
            current_batch.append(gn)
            current_tokens += ct_tokens

    if current_batch:
        batches.append(current_batch)

    return batches


def _anonymize_records(filtered_records: Dict[str, List[dict]]) -> Tuple[Dict[str, List[dict]], Dict[str, str]]:
    """Replace real group names with anonymous cluster IDs.

    Args:
        filtered_records: Dict of group_name -> list of DEG records

    Returns:
        Tuple of (anonymized filtered_records, anon_map: "CLUSTER N" -> real group_name)
    """
    anon_records = {}
    anon_map = {}  # "CLUSTER N" -> real group_name
    for i, gn in enumerate(sorted(filtered_records.keys()), 1):
        anon_id = f"CLUSTER {i}"
        anon_map[anon_id] = gn
        anon_records[anon_id] = filtered_records[gn]
    return anon_records, anon_map


def select_markers_from_deg(
    deg_records: List[dict],
    top_n: int,
    species: str,
    verbose: bool = False,
    metrics_path: Optional[Path] = None,
    command: str = "",
    anonymous: bool = False,
) -> List[dict]:
    """Select marker genes from DEG data using LLM.

    If the prompt exceeds ~160K tokens, auto-splits into multiple calls
    by batching cell types.

    Args:
        deg_records: List of evidence dicts from DEG extraction
        top_n: Number of top DEGs per cell type to show the LLM
        species: Species Latin name
        verbose: Whether to print progress
        metrics_path: Optional path to save LLM metrics
        command: CLI command string for metrics
        anonymous: Whether to anonymize cell type names as "CLUSTER N"

    Returns:
        List of evidence dicts with source_type="selected"
    """
    # Prepare DEG data
    _, filtered_records = prepare_deg_data_for_prompt(deg_records, top_n)

    if verbose:
        total_genes = sum(len(recs) for recs in filtered_records.values())
        print(f"\n   Prepared {len(filtered_records)} cell types, {total_genes} genes (top {top_n})")
        for gn, recs in sorted(filtered_records.items()):
            print(f"     - {gn}: {len(recs)} genes")

    # Anonymize if requested
    anon_map = None  # "CLUSTER N" -> real group_name
    if anonymous:
        filtered_records, anon_map = _anonymize_records(filtered_records)
        if verbose:
            print(f"   Anonymized {len(anon_map)} cell types as CLUSTER 1..{len(anon_map)}")

    # Estimate tokens and decide on batching
    all_cts = sorted(filtered_records.keys())
    total_est = _estimate_tokens(filtered_records, all_cts)

    if verbose:
        print(f"   Estimated input: ~{total_est:,} tokens")

    batches = _build_batches(filtered_records)

    if len(batches) > 1 and verbose:
        print(f"   Splitting into {len(batches)} batches to fit context window")

    # Run each batch
    prompt_name = "select_markers_anonymous" if anonymous else "select_markers"
    template = load_prompt_template(prompt_name)
    all_selections = []
    all_messages = []
    total_time = 0.0

    for batch_idx, batch_cts in enumerate(batches):
        batch_est = _estimate_tokens(filtered_records, batch_cts)
        if verbose and len(batches) > 1:
            print(f"\n   Batch {batch_idx + 1}/{len(batches)}: {len(batch_cts)} cell types, ~{batch_est:,} tokens")

        deg_text = _format_deg_sections(filtered_records, batch_cts)
        prompt = template.format(species=species, deg_data=deg_text)

        with Timer() as timer:
            result, message = call_claude_generate(
                prompt, SelectionResult, verbose=verbose
            )

        all_selections.extend(result.selections)
        all_messages.append(message)
        total_time += timer.elapsed

    # Save metrics from last message (or aggregated)
    if metrics_path and all_messages:
        last_msg = all_messages[-1]
        save_metrics(
            output_path=metrics_path,
            model=config.anthropic_model,
            message=last_msg,
            processing_time_sec=total_time,
            num_extractions=len(all_selections),
            command=command,
        )
        if verbose:
            print(f"   Saved metrics to: {metrics_path}")

    # Build lookup from filtered DEG data: (group_name_in_prompt, feature_name) -> record
    # When anonymous, keys use "CLUSTER N"; when named, keys use real group_name
    deg_lookup: Dict[Tuple[str, str], dict] = {}
    for gn, recs in filtered_records.items():
        for rec in recs:
            fn = (rec.get("feature_name") or "").strip().upper()
            deg_lookup[(gn, fn)] = rec

    # Validate selections and build evidence records
    anon_suffix = ":anon" if anonymous else ""
    source_id = f"selected:{config.anthropic_model}:top{top_n}{anon_suffix}"
    records = []
    discarded = 0

    for sel in all_selections:
        sel_gn = sel.group_name.strip().upper()
        fn = sel.feature_name.strip().upper()
        key = (sel_gn, fn)

        if key not in deg_lookup:
            discarded += 1
            if verbose:
                print(f"   Discarded (not in DEG top {top_n}): {sel_gn} / {fn}")
            continue

        deg_rec = deg_lookup[key]

        # Map back to real group name if anonymous
        real_gn = anon_map[sel_gn] if anon_map and sel_gn in anon_map else sel_gn

        records.append({
            "organism": species,
            "group_label": real_gn,
            "group_name": real_gn,
            "group_id": None,
            "feature_label": fn,
            "feature_name": fn,
            "feature_id": None,
            "source_type": "selected",
            "source_rationale": sel.rationale.strip(),
            "source_id": source_id,
            "data_id": deg_rec.get("data_id"),
            "metrics_pcorr": deg_rec.get("metrics_pcorr"),
            "metrics_logfc": deg_rec.get("metrics_logfc"),
            "metrics_rank": deg_rec.get("metrics_rank"),
        })

    if verbose:
        print(f"   Selected {len(records)} markers ({discarded} discarded)")

    return records
