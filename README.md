# mrkr

`mrkr` extracts cell type marker evidence from scientific manuscripts and anchors each term to a stable identifier.

The primary output is one `paper.onto.json` file per paper. This file is a machine-readable ledger between manuscript statements and ontology terms, similar to how [`span`](https://github.com/pachterlab/span) records the alignment between a mathematics paper and its Lean formalization.

## Workflow

The workflow has a hard boundary between language-model output and programmatic grounding.

```bash
# 1. Extract verbatim evidence and normalized labels. The LLM cannot emit identifiers.
uv run mrkr extract \
  --manuscript manuscript.md \
  --source-id doi:10.x/example \
  --output paper.claims.json

# 2. Assign Ensembl, Cell Ontology, and UBERON identifiers.
uv run mrkr ground paper.claims.json \
  --manuscript manuscript.md \
  --organism homo_sapiens \
  --output paper.onto.json

# 3. Verify source hashes, offsets, term cardinality, and grounding state.
uv run mrkr validate paper.onto.json \
  --manuscript manuscript.md
```

`extract` and `ground` write atomically. They do not leave an output file when validation fails.
`extract` records how many raw claims were retained, excluded because no gene was named, or
re-anchored to an exact manuscript quotation. It also records the model, response ID, and prompt
template digest used for the extraction.

## Claim document

Each claim contains one target cell type, at least one marker gene, and any comparison or tissue terms stated by the authors.

```json
{
  "schema_version": "mrkr.onto.v1",
  "source": {
    "id": "doi:10.x/example",
    "sha256": "sha256:..."
  },
  "grounding": {
    "genes": {
      "provider": "offline-gene-map",
      "organism": "homo_sapiens",
      "sha256": "sha256:..."
    },
    "ontology_service": {
      "provider": "OLS4",
      "endpoint": "https://www.ebi.ac.uk/ols4/api/v2/tag_text",
      "queries": [
        {
          "query": "conventional dendritic cell type 1",
          "ontology": "cl",
          "retrieved_at": "2026-07-21T00:00:00+00:00",
          "response_sha256": "sha256:...",
          "ontology_term": "CL:0000990",
          "exact": true
        },
        {
          "query": "myeloid cell",
          "ontology": "cl",
          "retrieved_at": "2026-07-21T00:00:00+00:00",
          "response_sha256": "sha256:...",
          "ontology_term": "CL:0000763",
          "exact": false
        }
      ]
    }
  },
  "claims": [
    {
      "claim_id": "claim:...",
      "span_literal": "cDC1 cells were distinguished from other myeloid cells by XCR1.",
      "span_offset": [120, 183],
      "summary": "conventional dendritic cell type 1 is distinguished from myeloid cell by XCR1.",
      "terms": [
        {
          "sub_span": "cDC1",
          "sub_offset": [0, 4],
          "normalized_label": "conventional dendritic cell type 1",
          "term_type": "celltype",
          "provenance": "explicit",
          "ontology_term": "CL:0000990",
          "exact": true
        },
        {
          "sub_span": "other myeloid cells",
          "sub_offset": [35, 54],
          "normalized_label": "myeloid cell",
          "term_type": "comparison",
          "provenance": "explicit",
          "ontology_term": "CL:0000763",
          "exact": false
        },
        {
          "sub_span": "XCR1",
          "sub_offset": [58, 62],
          "normalized_label": "XCR1",
          "term_type": "gene",
          "provenance": "explicit",
          "ontology_term": "ENSG00000173578",
          "exact": true,
          "direction": "positive"
        }
      ]
    }
  ]
}
```

Grounding has three explicit outcomes:

- `ontology_term` is set and `exact=true`: the full normalized label matched.
- `ontology_term` is set and `exact=false`: a broader ontology term matched.
- `ontology_term=null` and `exact=null`: the term remains unresolved.

An ontology service failure stops the command. It is not recorded as an unresolved term.

## Validation contract

A valid claim document satisfies these checks:

1. The source SHA-256 digest matches the manuscript.
2. Each evidence span is selected exactly by its manuscript offset.
3. Each explicit term span is selected exactly within the evidence span.
4. Each claim has exactly one target cell type and at least one marker gene.
5. Each marker gene is named explicitly in the evidence span.
6. Only genes have a positive or negative direction.
7. Grounded terms match recorded gene-map or ontology-query metadata.

`mrkr validate --report report.json` writes a machine-readable error report and exits with status 1 when any check fails.

## Installation

```bash
git clone https://github.com/sbooeshaghi/mrkr.git
cd mrkr
uv sync --locked
```

Create `.env` from `.env.example`, then set `ANTHROPIC_API_KEY`.

```bash
cp .env.example .env
```

Cell type, comparison, and tissue grounding uses the EBI Ontology Lookup Service. The packaged
gene map is human-only. For another organism, pass its name and a versioned two-column map with
`--organism` and `--gene-map`. A document must use one organism for gene grounding.
Each ontology request is recorded with its query, retrieval time, selected identifier, and response
digest. The resulting `paper.onto.json` is therefore the versioned result used downstream even if
the remote ontology service changes later.

## Development

```bash
uv sync --locked --extra dev
uv run --locked --extra dev pytest -q
uv run --locked --extra dev ruff check mrkr tests
```

The repository supports Python 3.10 through 3.13.

## License

MIT
