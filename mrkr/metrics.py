"""Metrics tracking for LLM usage and costs."""

import json
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Any

from anthropic.types import Message

try:
    from litellm import model_cost
    HAS_LITELLM = True
except ImportError:
    HAS_LITELLM = False


def calculate_anthropic_cost(message: Message, model_key: str) -> Dict[str, Any]:
    """Calculate cost breakdown for Anthropic API call.

    Args:
        message: Anthropic Message object with usage information
        model_key: Model identifier (e.g., 'claude-sonnet-4-5-20250929')

    Returns:
        Dictionary with token counts and costs
    """
    u = message.usage
    input_tokens = int(u.input_tokens)
    output_tokens = int(u.output_tokens)
    cache_read_tokens = int(getattr(u, "cache_read_input_tokens", 0))
    cache_creation_tokens = int(getattr(u, "cache_creation_input_tokens", 0))
    cached_tokens = cache_read_tokens + cache_creation_tokens
    total_tokens = input_tokens + output_tokens

    # Calculate costs if litellm is available
    input_cost = 0.0
    cached_cost = 0.0
    output_cost = 0.0
    total_cost = 0.0

    if HAS_LITELLM and model_key in model_cost:
        prices = model_cost[model_key]
        limit = prices.get("max_input_tokens", 200_000)

        # Split tokens into base (≤200k) and above (>200k)
        def split(n):
            return (min(n, limit), max(n - limit, 0))

        in_base, in_above = split(input_tokens)
        out_base, out_above = split(output_tokens)
        cr_base, cr_above = split(cache_read_tokens)
        cc_base, cc_above = split(cache_creation_tokens)

        # Prices
        in_p = prices["input_cost_per_token"]
        in_p_hi = prices.get("input_cost_per_token_above_200k_tokens", in_p)
        out_p = prices["output_cost_per_token"]
        out_p_hi = prices.get("output_cost_per_token_above_200k_tokens", out_p)
        cr_p = prices.get("cache_read_input_token_cost", 0.0)
        cr_p_hi = prices.get("cache_read_input_token_cost_above_200k_tokens", cr_p)
        cc_p = prices.get("cache_creation_input_token_cost", 0.0)
        cc_p_hi = prices.get("cache_creation_input_token_cost_above_200k_tokens", cc_p)

        # Costs
        input_cost = in_base * in_p + in_above * in_p_hi
        cached_cost = cr_base * cr_p + cr_above * cr_p_hi + cc_base * cc_p + cc_above * cc_p_hi
        output_cost = out_base * out_p + out_above * out_p_hi
        total_cost = input_cost + cached_cost + output_cost

    return {
        "input_tokens": input_tokens,
        "cached_tokens": cached_tokens,
        "output_tokens": output_tokens,
        "total_tokens": total_tokens,
        "input_cost": input_cost,
        "cached_cost": cached_cost,
        "output_cost": output_cost,
        "total_cost": total_cost,
    }


def save_metrics(
    output_path: Path,
    model: str,
    message: Message,
    processing_time_sec: float,
    num_extractions: int = 0,
) -> None:
    """Save metrics to JSON file.

    Args:
        output_path: Path to save metrics JSON file
        model: Model identifier string
        message: Anthropic Message object with usage information
        processing_time_sec: Time taken for the LLM call in seconds
        num_extractions: Number of marker extractions returned
    """
    # Calculate cost breakdown
    cost_breakdown = calculate_anthropic_cost(message, model)

    # Build full metrics dictionary
    metrics = {
        "model": model,
        "created_at": datetime.utcnow().isoformat() + 'Z',
        "timestamp": int(datetime.now().timestamp()),
        "num_extractions": num_extractions,
        "processing_time_sec": round(processing_time_sec, 2),
        **cost_breakdown,
    }

    # Write to file
    output_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")


class Timer:
    """Context manager for timing operations."""

    def __init__(self):
        self.start_time = None
        self.end_time = None
        self.elapsed = None

    def __enter__(self):
        self.start_time = time.time()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end_time = time.time()
        self.elapsed = self.end_time - self.start_time
        return False
