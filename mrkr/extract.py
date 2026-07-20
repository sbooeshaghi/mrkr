"""Main extraction orchestration for mrkr."""

import hashlib
from pathlib import Path
from typing import Dict, List, Optional

from .config import config
from .deg import extract_from_deg_files
from .llm import extract_claims_from_text_and_images, extract_from_text_and_images
from .metrics import save_metrics, Timer
from .verify import verify_extractions


def extract_markers(
    manuscript_path: Optional[Path],
    figure_paths: Optional[List[Path]],
    deg_paths: Optional[List[Path]],
    species: Optional[str],
    verbose: bool = False,
    metrics_path: Optional[Path] = None,
    verify: bool = True,
    known_cell_types: Optional[Dict[str, List[str]]] = None,
    command: str = "",
) -> List[dict]:
    """
    Main extraction orchestration.

    Workflows:
    1. mrkr --manuscript              → text extractions
    2. mrkr --manuscript --figures    → text+image extractions (joint call)
    3. mrkr --manuscript --deg        → text + DEG with cell type matching
    4. mrkr --manuscript --figures --deg → text+image (joint) + DEG with matching
    5. mrkr --deg                    → DEG only

    Args:
        manuscript_path: Path to manuscript markdown file
        figure_paths: List of figure image paths
        deg_paths: List of DEG table paths
        species: Species name (required for DEG-only)
        verbose: Whether to print progress
        metrics_path: Optional path to save LLM metrics JSON

    Returns:
        List of evidence dictionaries (uniform structure)
    """
    results = []
    deg_data_ids = {}  # Map group_name → list of data_ids
    manuscript_text = ""

    # STEP 1: Process DEG files if provided
    if deg_paths:
        if verbose:
            print(f"\n📊 Processing {len(deg_paths)} DEG file(s)...")

        # Extract from DEG
        deg_results = extract_from_deg_files(deg_paths, species, verbose)
        results.extend(deg_results)

        # Build mapping of cell types to data_ids for later matching
        # Also build per-source cell type dict for prompt context
        from collections import defaultdict
        deg_by_source = defaultdict(set)
        for record in deg_results:
            group_name = record["group_name"]
            data_id = record["data_id"]

            if group_name not in deg_data_ids:
                deg_data_ids[group_name] = []

            if data_id not in deg_data_ids[group_name]:
                deg_data_ids[group_name].append(data_id)

            deg_by_source[data_id].add(group_name)

        deg_ct_dict = {did: sorted(gns) for did, gns in sorted(deg_by_source.items())}

        # Merge DEG cell types with any provided known_cell_types
        if known_cell_types:
            for did, gns in deg_ct_dict.items():
                if did in known_cell_types:
                    known_cell_types[did] = sorted(set(known_cell_types[did]) | set(gns))
                else:
                    known_cell_types[did] = gns
        else:
            known_cell_types = deg_ct_dict

        if verbose:
            print(f"   ✓ Extracted {len(deg_results)} marker genes from DEG tables")
            print(f"   ✓ Found {len(known_cell_types)} unique cell types")

    # STEP 2: Process manuscript and/or figures (joint call if both provided)
    if manuscript_path or figure_paths:
        if verbose:
            if manuscript_path and figure_paths:
                print(
                    f"\n📄 Processing manuscript + {len(figure_paths)} figure(s) together..."
                )
            elif manuscript_path:
                print(f"\n📄 Processing manuscript...")
            else:
                print(f"\n🖼️  Processing {len(figure_paths)} figure(s)...")

        source_id = ""

        if manuscript_path:
            manuscript_text = manuscript_path.read_text(encoding="utf-8")
            source_id = manuscript_path.name

            if verbose:
                char_count = len(manuscript_text)
                print(f"   ✓ Read manuscript: {char_count:,} characters")

        if figure_paths and not manuscript_path:
            # Figures only (no manuscript)
            source_id = ", ".join(f.name for f in figure_paths)

        # Extract from text and/or images (joint call) with timing
        with Timer() as timer:
            text_image_results, message = extract_from_text_and_images(
                manuscript_text=manuscript_text,
                image_paths=figure_paths if figure_paths else None,
                source_id=source_id,
                known_cell_types=known_cell_types,
                verbose=verbose,
            )

        # Save metrics if requested
        if metrics_path:
            save_metrics(
                output_path=metrics_path,
                model=config.anthropic_model,
                message=message,
                processing_time_sec=timer.elapsed,
                num_extractions=len(text_image_results),
                command=command,
            )
            if verbose:
                print(f"   💾 Saved metrics to: {metrics_path}")

        # Post-process: fill in data_id from DEG mapping if LLM didn't set it
        if known_cell_types and deg_data_ids:
            for record in text_image_results:
                if record.get("data_id"):
                    continue  # LLM already assigned data_id
                group_name = record["group_name"]
                if group_name in deg_data_ids:
                    record["data_id"] = deg_data_ids[group_name][0]

        results.extend(text_image_results)

    # Validate we have some results
    if not results:
        raise ValueError(
            "No extractions found. Please check your inputs and try again."
        )

    if verbose:
        print(f"\n✅ Total: {len(results)} marker gene associations extracted")

    # STEP 3: Verify extractions against manuscript text
    if verify and manuscript_text:
        if verbose:
            print(f"\n🔍 Verifying extractions against manuscript...")
        results = verify_extractions(manuscript_text, results, verbose=verbose)

    return results


def extract_claims(
    manuscript_path: Optional[Path],
    figure_paths: Optional[List[Path]] = None,
    verbose: bool = False,
    metrics_path: Optional[Path] = None,
    validate: bool = True,
    command: str = "",
) -> List[dict]:
    """
    Extract marker CLAIM objects from a manuscript (± figures) — the current mrkr format.

    Returns a list of grounded-ready claim dicts: {paper_id, source_hash, span_literal,
    span_offset, summary, terms:[{sub_span, sub_offset, normalized_label, term_type,
    ontology_term(None), exact(None), direction, provenance}]}. Ids are filled by `mrkr ground`.
    """
    if not manuscript_path and not figure_paths:
        raise ValueError("extract_claims requires a manuscript and/or figures")

    manuscript_text = manuscript_path.read_text(encoding="utf-8") if manuscript_path else ""
    source_id = manuscript_path.name if manuscript_path else ", ".join(f.name for f in figure_paths)
    source_hash = ("sha256:" + hashlib.sha256(manuscript_text.encode("utf-8")).hexdigest()
                   if manuscript_text else None)

    if verbose:
        print(f"\n📄 Extracting claims from {source_id}...")

    with Timer() as timer:
        raw_claims, message = extract_claims_from_text_and_images(
            manuscript_text=manuscript_text,
            image_paths=figure_paths if figure_paths else None,
            source_id=source_id,
            verbose=verbose,
        )

    if metrics_path:
        save_metrics(output_path=metrics_path, model=config.anthropic_model, message=message,
                     processing_time_sec=timer.elapsed, num_extractions=len(raw_claims),
                     command=command)

    claims = []
    for c in raw_claims:
        span = c["span_literal"]
        span_offset = None
        if manuscript_text:
            i = manuscript_text.find(span)
            span_offset = [i, i + len(span)] if i >= 0 else None
        terms = []
        for t in c["terms"]:
            ss = t.get("sub_span")
            sub_offset = None
            if ss:
                j = span.find(ss)
                sub_offset = [j, j + len(ss)] if j >= 0 else None
            term = {
                "sub_span": ss,
                "sub_offset": sub_offset,
                "normalized_label": t["normalized_label"],
                "term_type": t["term_type"],
                "ontology_term": None,
                "exact": None,
                "provenance": "explicit" if ss else "implicit",
            }
            if t["term_type"] == "gene":
                term["direction"] = t.get("direction", "positive")
            terms.append(term)
        claims.append({
            "paper_id": source_id or None,
            "source_hash": source_hash,
            "span_literal": span,
            "span_offset": span_offset,
            "summary": c["summary"],
            "terms": terms,
        })

    if validate and manuscript_text:
        from .ground import validate as validate_claims
        r = validate_claims(claims, manuscript_text)
        if verbose:
            print(f"   🔍 validation — span {r['span_ok']}/{r['span_total']}, "
                  f"sub_span {r['sub_ok']}/{r['sub_total']}, label {r['label_ok']}/{r['label_total']}")

    if verbose:
        print(f"\n✅ Extracted {len(claims)} marker claims")

    return claims
