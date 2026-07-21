---
name: mrkr-authoring
description: Extract, validate, and improve source-grounded cell type marker statements in scientific manuscripts with mrkr. Use for manuscript authoring, marker-evidence audits, paper.onto.json generation, or checks that marker statements retain their target, comparison, tissue, direction, and verbatim source evidence.
---

# mrkr Authoring

Use `mrkr` to connect each reported marker to an exact manuscript span and typed biological terms. Treat `paper.onto.json` as the source of truth for the audit.

## Workflow

1. Locate the manuscript text and the `mrkr` project.
2. If `paper.onto.json` exists, validate it before reviewing the prose.
3. Otherwise, extract claims, ground terms, and validate the result.
4. Review the validated claims for missing biological context.
5. Edit manuscript prose only when the user requests edits.
6. Rerun extraction and validation after an edit that changes marker statements.

```bash
uv run --project /path/to/mrkr --locked mrkr extract \
  --manuscript manuscript.md \
  --source-id doi:10.x/example \
  --output paper.claims.json

uv run --project /path/to/mrkr --locked mrkr ground paper.claims.json \
  --manuscript manuscript.md \
  --organism homo_sapiens \
  --output paper.onto.json

uv run --project /path/to/mrkr --locked mrkr validate paper.onto.json \
  --manuscript manuscript.md \
  --report validation.json
```

Do not continue from an invalid document. Report the validation errors first.
Determine the paper organism from the manuscript or author metadata before grounding. Do not
guess when it is ambiguous. The packaged gene map supports only `homo_sapiens`; another organism
requires a versioned `--gene-map`.

## Authoring Review

For each claim, check these questions:

- Is the target cell population explicit?
- Are all reported positive and negative markers represented?
- Is the comparison population stated when the paper supports one?
- Is the relevant tissue, disease, perturbation, or assay context stated nearby?
- Does the wording distinguish a reported marker from an inferred association?
- Can a reader find the evidence from the stored span and offsets?

Treat a summary as incomplete when its evidence contains a qualifier that changes where or how the
marker applies but the summary omits it. Check comparison populations, species, tissue, disease or
perturbation, assay or validation method, direction, and quantitative restrictions such as a
reported fraction. Do not require background details that do not change the marker statement.

Classify findings as:

- `error`: the summary changes the target, marker, direction, or stated evidence.
- `missing qualifier`: the summary omits a supported condition that limits interpretation or reuse.
- `unresolved`: the source does not provide enough information to complete the statement.

Do not invent missing context. Record it as an authoring finding outside the onto file when the manuscript does not state it. Do not create an implicit term to represent absent information.

Report each actionable item with the claim ID, exact evidence span, missing or ambiguous field, and a concise proposed edit. Preserve the authors' biological meaning and terminology.

## Grounding Rules

- The language model emits spans and labels only.
- Programmatic grounding assigns Ensembl, Cell Ontology, and UBERON identifiers. The onto file records the grounding provider because remote ontology results can change.
- `exact=true` means that the full normalized label matched.
- `exact=false` means that only a broader ontology term matched.
- Null identifier and null exactness mean unresolved.
- Never treat a service failure as an unresolved term.

Read [references/claim-format.md](references/claim-format.md) when inspecting or modifying the JSON format.

## Cost Control

Rerun extraction after an accepted manuscript edit so the onto file remains synchronized. Before processing more than one manuscript, identify the corpus runner, perform its dry run, and obtain approval for the reported scope.

Never print API keys or include `.env` in generated artifacts.
