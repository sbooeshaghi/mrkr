"""LLM integration for marker gene extraction using Claude."""

import json
from pathlib import Path
from typing import List, Optional, Tuple

import anthropic
from anthropic.types import Message

from .config import config
from .image import process_image_for_api
from .models import ExtractionsResult


def load_prompt_template(template_name: str) -> str:
    """
    Load prompt template from prompts/ folder.

    Args:
        template_name: Name of template file (without .txt extension)

    Returns:
        Template content as string
    """
    prompts_dir = Path(__file__).parent / "prompts"
    template_path = prompts_dir / f"{template_name}.txt"

    if not template_path.exists():
        raise FileNotFoundError(f"Prompt template not found: {template_path}")

    return template_path.read_text(encoding="utf-8")


def extract_json_from_response(response_text: str) -> str:
    """Extract JSON from LLM response, handling markdown code fences.

    Args:
        response_text: Raw response text from LLM

    Returns:
        Cleaned JSON string
    """
    text = response_text.strip()

    # Check for markdown code fences
    if "```" in text:
        lines = text.split("\n")
        start_idx = None
        end_idx = None

        # Find the first opening fence
        for i, line in enumerate(lines):
            if line.strip().startswith("```"):
                start_idx = i + 1
                break

        # Find the closing fence
        if start_idx is not None:
            for i in range(start_idx, len(lines)):
                if lines[i].strip() == "```":
                    end_idx = i
                    break

            if end_idx is not None and start_idx < end_idx:
                # Extract content between fences
                json_lines = lines[start_idx:end_idx]
                return "\n".join(json_lines).strip()

    # If no code fences or extraction failed, look for JSON object boundaries
    # Find the first '{' and last '}'
    first_brace = text.find("{")
    last_brace = text.rfind("}")

    if first_brace != -1 and last_brace != -1 and first_brace < last_brace:
        return text[first_brace:last_brace + 1]

    # Return as-is if no clear JSON structure found
    return text


def call_claude_json(
    prompt: str,
    image_paths: Optional[List[Path]] = None,
    verbose: bool = False,
) -> Tuple[ExtractionsResult, Message]:
    """
    Call Claude API with prompt and optional images, expecting JSON response.

    Args:
        prompt: The prompt text
        image_paths: Optional list of image file paths for vision mode
        verbose: Whether to print token usage

    Returns:
        Tuple of (ExtractionsResult, Message) for metrics tracking
    """
    client = anthropic.Anthropic(api_key=config.anthropic_api_key)

    # Prepare message content
    content_parts = []

    # Add images first if provided (resize and encode for API)
    if image_paths:
        for img_path in image_paths:
            # Process image (resize to Claude limits and encode as base64)
            result = process_image_for_api(img_path, max_size=1568, verbose=verbose)

            if result is None:
                if verbose:
                    print(f"   ⚠️  Skipping image (processing failed): {img_path.name}")
                continue

            image_data, media_type = result

            content_parts.append({
                "type": "image",
                "source": {
                    "type": "base64",
                    "media_type": media_type,
                    "data": image_data,
                },
            })

    # Add text prompt
    content_parts.append({"type": "text", "text": prompt})

    # Call Claude
    if verbose:
        img_count = len(image_paths) if image_paths else 0
        print(f"   🤖 Calling Claude ({img_count} images)...")

    message = client.messages.create(
        model=config.anthropic_model,
        max_tokens=65536,  # Maximum for Claude Sonnet 4.5 (64k output tokens)
        temperature=0.0,  # Deterministic output for consistent JSON structure
        messages=[{"role": "user", "content": content_parts}],
    )

    response_text = message.content[0].text

    if verbose:
        print(
            f"   ✓ Response: {message.usage.input_tokens} tokens in, "
            f"{message.usage.output_tokens} tokens out"
        )

    # Parse JSON response
    try:
        # Extract JSON from response (handle markdown code fences)
        json_text = extract_json_from_response(response_text)

        if not json_text or not json_text.strip():
            raise ValueError(
                f"Failed to extract JSON from response. Response text: {response_text[:1000]}"
            )

        result = ExtractionsResult.model_validate_json(json_text)
        return result, message

    except Exception as e:
        if verbose:
            print(f"   ❌ Failed to parse JSON response")
            print(f"   Raw response:\n{response_text[:500]}...")
        raise ValueError(f"Failed to parse Claude's response as JSON: {e}")


def extract_from_text_and_images(
    manuscript_text: str,
    image_paths: Optional[List[Path]],
    source_id: str,
    known_cell_types: Optional[List[str]] = None,
    verbose: bool = False,
) -> Tuple[List[dict], Message]:
    """
    Extract marker genes from manuscript text and/or images using Claude.

    If both text and images are provided, they are processed together in a single call.

    Args:
        manuscript_text: Manuscript content (can be empty if only images)
        image_paths: List of figure image paths (optional)
        source_id: Identifier for the source (e.g., filename)
        known_cell_types: Optional list of known cell type names from DEG
        verbose: Whether to print progress

    Returns:
        Tuple of (extractions list, Message) for metrics tracking
    """
    # Choose prompt template based on inputs
    if image_paths and known_cell_types:
        template_name = "extract_text_and_images_with_deg"
    elif image_paths:
        template_name = "extract_text_and_images"
    elif known_cell_types:
        template_name = "extract_text_with_deg"
    else:
        template_name = "extract_text"

    # Load and format prompt
    prompt_template = load_prompt_template(template_name)

    # Prepare template variables
    template_vars = {"manuscript_text": manuscript_text}

    # Add known cell types if provided
    if known_cell_types:
        cell_types_list = "\n".join(f"  - {ct}" for ct in sorted(known_cell_types))
        template_vars["known_cell_types"] = cell_types_list

    # Inject variables into template
    prompt = prompt_template.format(**template_vars)

    # Call Claude (with images if provided)
    result, message = call_claude_json(prompt, image_paths=image_paths, verbose=verbose)

    # Convert to dictionaries (match Evidence structure)
    extractions = []
    for ext in result.extractions:
        extractions.append({
            "organism": ext.organism.strip(),
            "group_label": ext.group_label.strip(),
            "group_name": ext.group_name.strip().upper(),
            "group_id": None,
            "feature_label": ext.feature_label.strip(),
            "feature_name": ext.feature_name.strip().upper(),
            "feature_id": None,
            "source_type": ext.source_type,
            "source_rationale": ext.source_rationale.strip(),
            "source_id": source_id,
            "data_id": None,  # Will be set during post-processing if DEG provided
            "metrics_pcorr": None,
            "metrics_logfc": None,
            "metrics_rank": None,
        })

    if verbose:
        print(f"   ✓ Extracted {len(extractions)} marker gene associations")

    return extractions, message
