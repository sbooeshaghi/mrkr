"""Anthropic integration for marker claim extraction."""

from __future__ import annotations

from pathlib import Path
from typing import TypeVar

import anthropic
from anthropic.types import Message
from pydantic import BaseModel

from .config import config
from .models import ClaimsResult

ResultModel = TypeVar("ResultModel", bound=BaseModel)


def load_prompt_template(template_name: str) -> str:
    """Load one packaged prompt template."""

    path = Path(__file__).parent / "prompts" / f"{template_name}.txt"
    if not path.is_file():
        raise FileNotFoundError(f"prompt template not found: {path}")
    return path.read_text(encoding="utf-8")


def extract_json_from_response(response_text: str) -> str:
    """Extract a JSON object from a plain or fenced model response."""

    text = response_text.strip()
    if "```" in text:
        lines = text.splitlines()
        opening = next(
            (index for index, line in enumerate(lines) if line.strip().startswith("```")),
            None,
        )
        if opening is not None:
            closing = next(
                (
                    index
                    for index in range(opening + 1, len(lines))
                    if lines[index].strip() == "```"
                ),
                None,
            )
            if closing is not None:
                return "\n".join(lines[opening + 1 : closing]).strip()

    first = text.find("{")
    last = text.rfind("}")
    return text[first : last + 1] if 0 <= first < last else text


def call_claude_json(
    prompt: str,
    result_model: type[ResultModel],
    verbose: bool = False,
) -> tuple[ResultModel, Message]:
    """Call Anthropic once and parse one typed JSON result."""

    client = anthropic.Anthropic(
        api_key=config.anthropic_api_key,
        timeout=config.timeout,
    )

    def stream(**extra) -> Message:
        with client.messages.stream(
            model=config.anthropic_model,
            max_tokens=32_000,
            messages=[{"role": "user", "content": prompt}],
            **extra,
        ) as response:
            return response.get_final_message()

    if verbose:
        print(f"Calling {config.anthropic_model}")
    try:
        message = stream(temperature=0.0)
    except anthropic.BadRequestError as error:
        if "temperature" not in str(error).lower():
            raise
        message = stream()

    response_text = next(
        (
            block.text
            for block in message.content
            if getattr(block, "type", None) == "text"
        ),
        "",
    )
    json_text = extract_json_from_response(response_text)
    if not json_text:
        raise ValueError("model response did not contain JSON")
    try:
        result = result_model.model_validate_json(json_text)
    except Exception:
        from json_repair import repair_json

        result = result_model.model_validate(repair_json(json_text, return_objects=True))

    if verbose:
        print(
            f"Received {message.usage.input_tokens} input tokens and "
            f"{message.usage.output_tokens} output tokens"
        )
    return result, message


def extract_claims_from_text(
    manuscript_text: str,
    organism_label: str,
    verbose: bool = False,
) -> tuple[list[dict], Message]:
    """Extract raw marker claim objects from complete manuscript text."""

    prompt = load_prompt_template("extract_claims").format(
        manuscript_text=manuscript_text,
        organism_label=organism_label,
    )
    result, message = call_claude_json(
        prompt,
        result_model=ClaimsResult,
        verbose=verbose,
    )

    claims: list[dict] = []
    for claim in result.claims:
        terms: list[dict] = []
        for extracted in claim.terms:
            term = {
                "sub_span": extracted.sub_span.strip()
                if extracted.sub_span
                else None,
                "normalized_label": extracted.normalized_label.strip(),
                "term_type": extracted.term_type,
            }
            if extracted.term_type == "gene":
                term["direction"] = extracted.direction
            terms.append(term)
        claims.append(
            {
                "span_literal": claim.span_literal,
                "summary": claim.summary.strip(),
                "terms": terms,
            }
        )
    return claims, message
