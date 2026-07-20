# mrkr - Cell Type Marker Gene Extraction

**mrkr** (marker) is a CLI tool that extracts cell type marker genes from scientific manuscripts and anchors them in ontologies.

## The flow: extract, then ground

The core is two steps with a hard line between them. The LLM proposes in natural language; the machine assigns the identifiers.

```bash
# 1. extract: LLM reads the paper -> claim objects (spans + normalized labels, NO ontology ids)
mrkr extract -m manuscript.md -o claims.json

# 2. ground: deterministic -> assign ids (gene -> Ensembl, celltype/comparison -> CL, tissue -> UBERON)
mrkr ground claims.json -o grounded.json -m manuscript.md
```

A claim is a verbatim sentence plus a flat list of grounded terms:

```json
{
  "span_literal": "cDC1 cells were distinguished from other myeloid cells by XCR1 and CLEC9A.",
  "summary": "Conventional dendritic cell type 1 is distinguished from other myeloid cells by XCR1 and CLEC9A.",
  "terms": [
    {"sub_span": "cDC1", "normalized_label": "conventional dendritic cell type 1", "term_type": "celltype", "ontology_term": "CL:0000990", "exact": true},
    {"sub_span": "XCR1", "normalized_label": "XCR1", "term_type": "gene", "ontology_term": "ENSG00000173578", "exact": true, "direction": "positive"}
  ]
}
```

The LLM never emits an id. Grounding assigns them: genes via the offline `gmap.txt`; cell types and tissue via OLS `tag_text` (singularized query, longest span-coverage → `exact`, coarse parents flagged `exact=false`). Every id traces to a verbatim span. See the LLMarkers design note (`llmarkers/docs/notes/mrkr_format_reframe_design_note.md`) for the full rationale.

## Features

- 📄 **Extract from manuscripts**: Parse markdown/text manuscripts to find cell type-marker gene associations
- 🖼️ **Extract from figures**: Analyze heatmaps, UMAPs, violin plots, and other figure types
- 📊 **Process DEG tables**: Handle Excel, CSV, and TSV differential expression tables with flexible column matching
- 🤖 **Claude-powered**: Uses Anthropic's Claude for intelligent extraction with species inference
- 🔄 **Unified output**: Consistent JSON format regardless of source type
- 🎯 **Smart matching**: When DEG tables are provided, matches manuscript mentions to DEG cell type names
- 🔎 **Grounded query**: Resolve natural-language marker-profile questions against an LLMarkers SQLite database

## Installation

```bash
# Using pip
pip install mrkr

# Using uv
uv pip install mrkr

# For development
git clone <repository>
cd mrid_new
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
uv pip install -e .
```

## Configuration

Create a `.env` file in your working directory or set environment variables:

```bash
ANTHROPIC_API_KEY=your_api_key_here
ANTHROPIC_MODEL=claude-sonnet-4-5-20250929
LLM_TIMEOUT=600.0
```

Get your API key from: https://console.anthropic.com/

## Usage

### Basic Examples

```bash
# Extract from manuscript only
mrkr -m manuscript.md -o markers.json

# Extract from manuscript + figures (processed together)
mrkr -m manuscript.md -f fig1.png -f fig2.png -o markers.json

# Extract from DEG table only
mrkr -d deg_table.xlsx -s homo_sapiens -o markers.json

# Full pipeline: manuscript + figures + DEG
mrkr -m manuscript.md -f fig1.png -f fig2.png -d deg.xlsx -o markers.json

# Multiple DEG files
mrkr -d deg1.xlsx -d deg2.csv -s homo_sapiens -o markers.json

# Verbose mode (show token usage and progress)
mrkr -m manuscript.md -d deg.xlsx -o markers.json -v
```

### Query LLMarkers

```bash
# Resolve a marker-pattern query against an LLMarkers SQLite database
mrkr query "TREM2+ macrophages in tumors" --db docs/llmarkers.sqlite

# Add explicit marker genes when the natural-language query does not name them all
mrkr query "exhausted T cells in melanoma" \
  --gene PDCD1 --gene HAVCR2 --gene LAG3 \
  --db docs/llmarkers.sqlite

# Save the full evidence card JSON
mrkr query "CCR8+ Tregs in tumors" --db docs/llmarkers.sqlite -o ccr8_treg_query.json
```

`mrkr query` uses an LLM by default to parse the natural-language prompt into a
structured query containing a cell type label, optional context, and marker genes
mapped to Ensembl IDs. Use `--parser heuristic` for an offline deterministic
parser. After parsing, retrieval is deterministic: marker matches use Ensembl
gene overlap, label matches use exact/partial token matching, and context
matches use the source and paper text stored in the database. The command returns
ranked evidence cards with profile IDs, paper titles, DOIs, shared genes, label
relations, source sentences, and conservative relationship calls.

### Options

```
-m, --manuscript PATH   Manuscript markdown/text file
-f, --figures PATH      Figure images (can specify multiple)
-d, --deg PATH          DEG table files (can specify multiple)
-o, --output PATH       Output JSON file (required)
-s, --species TEXT      Species Latin name (e.g., homo_sapiens)
-v, --verbose           Show detailed progress and token usage
--version               Show version
--help                  Show help message
```

## Input Formats

### Manuscripts
- **Formats**: `.md`, `.txt`
- **Content**: Scientific manuscripts in markdown or plain text
- **Processing**: Claude extracts cell type-marker gene mentions and infers species per extraction

### Figures
- **Formats**: `.png`, `.jpg`, `.jpeg`
- **Types**: Heatmaps, UMAPs, violin plots, dot plots, tables
- **Processing**: Claude vision analyzes figures for marker genes
- **Joint processing**: When manuscript + figures provided together, processed in single call

### DEG Tables
- **Formats**: `.xlsx`, `.csv`, `.tsv`
- **Multi-sheet**: Excel files process each sheet separately
- **Flexible columns**: Automatically detects column names (case-insensitive)

#### Required Columns
- **Cluster/Cell type**: `cluster`, `celltype`, `cell_type`, `cell-type`, `group`, `orig.ident`
- **Gene**: `gene`, `gene_name`, `genename`, `gene name`
- **Log fold change**: `logfc`, `log_fc`, `log2fc`, `avg_log2fc`, `avg_logFC`, etc.

#### Optional Columns
- **Adjusted p-value**: `p_corr`, `p_val_adj`, `pval_adj`, `padj`, etc.

## Output Format

All extractions use a uniform JSON structure:

```json
[
  {
    "organism": "homo_sapiens",
    "group_label": "naive CD4+ T cells",
    "group_name": "NAIVE CD4+ T CELL",
    "group_id": null,
    "feature_label": "CD4",
    "feature_name": "CD4",
    "feature_id": null,
    "source_type": "text",
    "source_rationale": "Clustering revealed naive CD4+ T cells expressing CD4 marker.",
    "source_id": "manuscript.md",
    "data_id": "deg_table#Sheet1",
    "metrics_pcorr": null,
    "metrics_logfc": null,
    "metrics_rank": null
  },
  {
    "organism": "homo_sapiens",
    "group_label": "B cell",
    "group_name": "B CELL",
    "group_id": null,
    "feature_label": "CD19",
    "feature_name": "CD19",
    "feature_id": null,
    "source_type": "deg",
    "source_rationale": "unfiltered",
    "source_id": "deg_table#Sheet1",
    "data_id": "deg_table#Sheet1",
    "metrics_pcorr": 0.0,
    "metrics_logfc": 2.5,
    "metrics_rank": 1
  }
]
```

### Field Descriptions

- **organism**: Species Latin name (e.g., `homo_sapiens`, `mus_musculus`)
- **group_label**: Exact cell type text as appears in source
- **group_name**: Normalized cell type name (UPPERCASE)
- **group_id**: Cell ontology ID (populated later, initially `null`)
- **feature_label**: Exact gene name as appears in source
- **feature_name**: Normalized gene name (UPPERCASE)
- **feature_id**: Gene ID like ENSEMBL (populated later, initially `null`)
- **source_type**: `"text"`, `"image"`, or `"deg"`
- **source_rationale**: Evidence text snippet or description
- **source_id**: Source filename
- **data_id**: DEG table identifier (for matching, `null` otherwise)
- **metrics_pcorr**: Adjusted p-value (DEG only, `null` otherwise)
- **metrics_logfc**: Log fold change (DEG only, `null` otherwise)
- **metrics_rank**: Rank within cell type by logFC (DEG only, `null` otherwise)

## Workflows

### 1. Manuscript Only
```bash
mrkr -m manuscript.md -o markers.json
```
- Extracts cell type-marker gene pairs from text
- Claude infers species for each extraction
- Normalizes cell type and gene names

### 2. Manuscript + Figures
```bash
mrkr -m manuscript.md -f fig1.png -f fig2.png -o markers.json
```
- Processes manuscript and figures **together** in single Claude call
- Extracts from both text and images
- `source_type` indicates whether from `"text"` or `"image"`

### 3. Manuscript + DEG
```bash
mrkr -m manuscript.md -d deg_table.xlsx -o markers.json
```
- Extracts cell types from DEG table
- Uses DEG cell type names as reference for manuscript extraction
- Claude matches manuscript mentions to DEG cell types
- Sets `data_id` field to link manuscript extractions to DEG tables
- Includes DEG extractions with metrics

### 4. DEG Only
```bash
mrkr -d deg_table.xlsx -s homo_sapiens -o markers.json
```
- Direct extraction from DEG tables
- Requires `--species` flag
- Computes ranks within cell types

## Species Support

mrkr uses Latin names with underscores for species:
- `homo_sapiens` (human)
- `mus_musculus` (mouse)
- `rattus_norvegicus` (rat)
- etc.

When processing manuscripts, Claude infers the species for each extraction independently, allowing for manuscripts that discuss multiple species.

## Error Handling

mrkr handles various edge cases:
- Missing or malformed columns in DEG tables
- Multiple Excel sheets (processed separately)
- Header row detection (tries up to 11 skip rows)
- Unicode in manuscripts
- Large manuscripts (handled by Claude's context window)
- Missing species (inferred from manuscript or required for DEG-only)

## Development

```bash
# Clone repository
git clone <repository>
cd mrid_new

# Create virtual environment
uv venv
source .venv/bin/activate

# Install in editable mode
uv pip install -e .

# Run tests (TODO)
pytest tests/

# Format code
black mrkr/
```

## Project Structure

```
mrkr/
├── mrkr/
│   ├── __init__.py          # Package initialization
│   ├── cli.py               # Click CLI interface
│   ├── config.py            # Configuration management
│   ├── models.py            # Pydantic models
│   ├── extract.py           # Main orchestration
│   ├── llm.py              # Claude API integration
│   ├── deg.py              # DEG file processing
│   ├── utils.py            # Utility functions
│   └── prompts/            # Prompt templates
│       ├── extract_text.txt
│       ├── extract_text_with_deg.txt
│       ├── extract_text_and_images.txt
│       └── extract_text_and_images_with_deg.txt
├── examples/                # Example files
├── tests/                   # Test suite
├── pyproject.toml          # Package configuration
├── README.md               # This file
└── .env.example            # Example environment config
```

## License

MIT

## Citation

If you use mrkr in your research, please cite:

```bibtex
@software{mrkr2024,
  title={mrkr: Cell Type Marker Gene Extraction Tool},
  author={mrkr developers},
  year={2024},
  url={https://github.com/...}
}
```

## Support

- Issues: https://github.com/.../issues
- Documentation: https://github.com/.../wiki
