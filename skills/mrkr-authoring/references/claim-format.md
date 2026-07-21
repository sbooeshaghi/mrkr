# Claim Format

`mrkr.claims.v1` is the ungrounded extraction document. `mrkr.onto.v1` is the grounded document. Both use one source envelope and a list of claims.

Each claim contains:

- `claim_id`: deterministic claim identifier.
- `span_literal` and `span_offset`: exact evidence in the source manuscript.
- `summary`: human-readable canonical rewrite used for review.
- `terms`: one target cell type, one or more genes, and optional comparison or tissue terms.

Each term contains:

- `sub_span` and `sub_offset`: exact surface text, or null for an implicit term.
- `normalized_label`: canonical text used for ontology grounding and comparison.
- `term_type`: `gene`, `celltype`, `comparison`, or `tissue`.
- `provenance`: `explicit` or `implicit`.
- `ontology_term`: CURIE or Ensembl identifier in an onto document; null if unresolved.
- `exact`: true for a full match, false for a coarse match, and null if unresolved.
- `direction`: `positive` or `negative` for genes only.

A valid claim has exactly one target cell type and at least one marker gene. Every marker gene must
be explicit in the evidence span. Comparisons are optional because many papers do not report them.
Their absence is a reporting limitation, not a schema error.

The source SHA-256 digest, evidence offsets, term offsets, cardinality, normalized labels,
direction, organism, and grounding state must all validate before an onto file is used downstream.
