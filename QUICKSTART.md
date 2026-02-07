# mrkr Quick Start Guide

## Installation

```bash
# Create a fresh virtual environment
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install click anthropic pydantic pandas openpyxl

# Install mrkr in development mode
pip install -e .
```

## Configuration

Create a `.env` file:

```bash
ANTHROPIC_API_KEY=sk-ant-...your-key-here...
```

Get your API key from: https://console.anthropic.com/

## Quick Test

```bash
# Check installation
mrkr --version
mrkr --help

# Test with a simple file (create test.md first)
echo "B cells express CD19 marker." > test.md
mrkr -m test.md -o output.json -v
```

## Common Commands

```bash
# Manuscript only
mrkr -m manuscript.md -o markers.json

# Manuscript + figures (joint processing)
mrkr -m manuscript.md -f fig1.png -f fig2.png -o markers.json

# DEG table only (requires species)
mrkr -d deg_table.xlsx -s homo_sapiens -o markers.json

# Full pipeline
mrkr -m manuscript.md -f fig1.png -d deg.xlsx -o markers.json

# Verbose mode
mrkr -m manuscript.md -o markers.json -v
```

## Expected Output Structure

```json
[
  {
    "organism": "homo_sapiens",
    "group_label": "B cells",
    "group_name": "B CELL",
    "group_id": null,
    "feature_label": "CD19",
    "feature_name": "CD19",
    "feature_id": null,
    "source_type": "text",
    "source_rationale": "B cells express CD19 marker.",
    "source_id": "test.md",
    "data_id": null,
    "metrics_pcorr": null,
    "metrics_logfc": null,
    "metrics_rank": null
  }
]
```

## Troubleshooting

### "ANTHROPIC_API_KEY required"
- Make sure `.env` file exists in your working directory
- Or set environment variable: `export ANTHROPIC_API_KEY=your-key`

### "No such file or directory"
- Use absolute paths or run from the correct directory
- Check file extensions (.md, .txt for manuscripts)

### "Could not find required columns"
- DEG tables must have: cluster/celltype, gene, logfc columns
- Check column names (case-insensitive matching)
- Try opening the file in Excel to verify format

### Import errors
- Make sure all dependencies are installed
- Run: `pip install click anthropic pydantic pandas openpyxl`

## Project Structure

```
mrkr/
├── mrkr/                    # Source code
│   ├── cli.py              # Command-line interface
│   ├── config.py           # Configuration (.env loading)
│   ├── models.py           # Data models
│   ├── extract.py          # Main orchestration
│   ├── llm.py             # Claude API calls
│   ├── deg.py             # DEG file processing
│   ├── utils.py           # Utilities
│   └── prompts/           # Prompt templates (4 .txt files)
├── examples/               # Example data
├── pyproject.toml         # Package config
├── README.md              # Full documentation
└── QUICKSTART.md          # This file
```

## Next Steps

1. Read the full [README.md](README.md) for detailed documentation
2. Try with your own manuscripts and DEG tables
3. Explore the prompt templates in `mrkr/prompts/`
4. Adjust species names as needed for your data
5. Post-process the output JSON as needed (add gene IDs, cell ontology IDs, etc.)

## Support

For issues or questions, please check:
- README.md for full documentation
- Prompt templates in mrkr/prompts/ to understand extraction logic
- Example outputs to verify format
