"""Extraction orchestration for one manuscript."""

from __future__ import annotations

import hashlib
from pathlib import Path

from .claims import assert_valid_document, make_claim_document, prepare_raw_claims
from .config import config
from .llm import extract_claims_from_text, load_prompt_template
from .metrics import Timer, save_metrics


def extract_claims(
    manuscript_path: Path,
    source_id: str | None = None,
    verbose: bool = False,
    metrics_path: Path | None = None,
    validate: bool = True,
    command: str = "",
) -> dict:
    """Extract and validate source-grounded marker claims from one manuscript."""

    manuscript_text = manuscript_path.read_text(encoding="utf-8")
    source_id = source_id or manuscript_path.name
    if verbose:
        print(f"Extracting marker claims from {source_id}")

    with Timer() as timer:
        raw_claims, message = extract_claims_from_text(
            manuscript_text=manuscript_text,
            verbose=verbose,
        )

    prepared_claims, preparation = prepare_raw_claims(manuscript_text, raw_claims)
    document = make_claim_document(
        source_id=source_id,
        manuscript_text=manuscript_text,
        raw_claims=prepared_claims,
    )
    prompt_template = load_prompt_template("extract_claims")
    document["extraction"] = {
        "provider": "anthropic",
        "model": config.anthropic_model,
        "response_id": getattr(message, "id", None),
        "prompt_template_sha256": "sha256:"
        + hashlib.sha256(prompt_template.encode("utf-8")).hexdigest(),
        "preparation": preparation,
    }

    if metrics_path:
        save_metrics(
            output_path=metrics_path,
            model=config.anthropic_model,
            message=message,
            processing_time_sec=timer.elapsed,
            num_extractions=len(document["claims"]),
            command=command,
        )

    report = assert_valid_document(document, manuscript_text) if validate else None
    if verbose and report is not None:
        print(
            f"Validated {report['n_claims']} claims; "
            f"re-anchored {preparation['reanchored_spans']}; "
            f"excluded {len(preparation['excluded_claims'])}"
        )
    return document
