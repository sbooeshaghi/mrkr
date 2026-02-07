"""Configuration management for mrkr."""

import os
from pathlib import Path


def load_env_file():
    """Load .env file from current directory or parent directories."""
    current = Path.cwd()
    for _ in range(4):
        env_file = current / ".env"
        if env_file.exists():
            with open(env_file) as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#") and "=" in line:
                        key, value = line.split("=", 1)
                        key = key.strip()
                        value = value.strip()
                        # Only set if not already in environment
                        if key not in os.environ:
                            os.environ[key] = value
            return

        parent = current.parent
        if parent == current:  # Reached root
            break
        current = parent


# Load .env file on import
load_env_file()


class Config:
    """Configuration for mrkr tool."""

    def __init__(self):
        """Initialize configuration from environment variables."""
        self.anthropic_api_key = os.getenv("ANTHROPIC_API_KEY", "")
        self.anthropic_model = os.getenv("ANTHROPIC_MODEL", "claude-sonnet-4-5-20250929")
        self.timeout = float(os.getenv("LLM_TIMEOUT", "600.0"))

    def validate(self):
        """Validate configuration.

        Raises:
            ValueError: If required configuration is missing
        """
        if not self.anthropic_api_key:
            raise ValueError(
                "ANTHROPIC_API_KEY environment variable is required. "
                "Set it with: export ANTHROPIC_API_KEY=your_api_key"
            )


# Global configuration instance
config = Config()
