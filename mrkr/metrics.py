"""Metrics tracking for language-model usage."""

import json
import time
from datetime import datetime
from pathlib import Path

from anthropic.types import Message


def usage_metrics(message: Message) -> dict[str, int]:
    """Return provider-reported token counts without embedding price data."""

    u = message.usage
    input_tokens = int(u.input_tokens)
    output_tokens = int(u.output_tokens)
    cache_read_tokens = int(getattr(u, "cache_read_input_tokens", 0))
    cache_creation_tokens = int(getattr(u, "cache_creation_input_tokens", 0))
    cached_tokens = cache_read_tokens + cache_creation_tokens

    return {
        "input_tokens": input_tokens,
        "cached_tokens": cached_tokens,
        "output_tokens": output_tokens,
        "total_tokens": input_tokens + output_tokens,
    }


def save_metrics(
    output_path: Path,
    model: str,
    message: Message,
    processing_time_sec: float,
    num_extractions: int = 0,
    command: str = "",
) -> None:
    """Save metrics to JSON file.

    Args:
        output_path: Path to save metrics JSON file
        model: Model identifier string
        message: Anthropic Message object with usage information
        processing_time_sec: Time taken for the LLM call in seconds
        num_extractions: Number of marker extractions returned
        command: Full CLI command string for reproducibility
    """
    metrics = {
        "command": command,
        "model": model,
        "created_at": int(datetime.now().timestamp()),
        "num_extractions": num_extractions,
        "processing_time_sec": round(processing_time_sec, 2),
        **usage_metrics(message),
    }

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
