"""Canonical claim documents and strict provenance validation."""

from __future__ import annotations

import copy
import hashlib
import re
from dataclasses import dataclass
from difflib import SequenceMatcher
from typing import Any, Iterable

from . import __version__

CLAIMS_SCHEMA = "mrkr.claims.v1"
ONTO_SCHEMA = "mrkr.onto.v1"
TERM_TYPES = {"gene", "celltype", "comparison", "tissue", "organism"}
GENE_DIRECTIONS = {"positive", "negative"}


@dataclass(frozen=True)
class ValidationIssue:
    """One machine-readable validation failure."""

    code: str
    path: str
    message: str

    def as_dict(self) -> dict[str, str]:
        return {"code": self.code, "path": self.path, "message": self.message}


class ClaimValidationError(ValueError):
    """Raised when a claim document violates the mrkr contract."""

    def __init__(self, report: dict[str, Any]):
        self.report = report
        errors = report.get("errors", [])
        preview = "; ".join(
            f"{item['path']}: {item['message']}" for item in errors[:3]
        )
        more = f"; {len(errors) - 3} more" if len(errors) > 3 else ""
        super().__init__(f"claim validation failed ({len(errors)} errors): {preview}{more}")


def sha256_text(text: str) -> str:
    """Return the canonical source digest used by mrkr documents."""

    return "sha256:" + hashlib.sha256(text.encode("utf-8")).hexdigest()


def _find_offset(container: str, value: str | None) -> list[int] | None:
    if not value:
        return None
    start = container.find(value)
    return [start, start + len(value)] if start >= 0 else None


def _all_offsets(container: str, value: str) -> list[int]:
    offsets: list[int] = []
    start = 0
    while True:
        index = container.find(value, start)
        if index < 0:
            return offsets
        offsets.append(index)
        start = index + 1


def _reanchor_span(manuscript_text: str, span: str) -> str | None:
    """Recover an exact source span when the model omitted text in its quotation."""

    if not span:
        return None
    if span in manuscript_text:
        return span

    # Require exact anchors at both ends. This handles omitted citations or
    # parentheticals without accepting a paraphrase as evidence.
    for width in (64, 48, 32, 24, 16, 12):
        if len(span) < width * 2:
            continue
        prefix = span[:width]
        suffix = span[-width:]
        starts = _all_offsets(manuscript_text, prefix)
        ends = _all_offsets(manuscript_text, suffix)
        candidates: list[str] = []
        max_length = max(len(span) * 3, len(span) + 2_000)
        for start in starts:
            for suffix_start in ends:
                end = suffix_start + len(suffix)
                if start < end and end - start <= max_length:
                    candidate = manuscript_text[start:end]
                    if SequenceMatcher(None, span, candidate).ratio() >= 0.80:
                        candidates.append(candidate)
        unique = list(dict.fromkeys(candidates))
        if len(unique) == 1:
            return unique[0]
    return None


def _term_position(span: str, sub_span: str, target_position: int | None) -> int:
    positions = _all_offsets(span, sub_span)
    if not positions:
        return len(span)
    if target_position is None:
        return positions[0]
    return min(positions, key=lambda position: abs(position - target_position))


def _reconstruct_span_from_terms(
    manuscript_text: str, claim: dict[str, Any]
) -> str | None:
    """Find the smallest ordered source interval containing explicit term anchors."""

    span = claim.get("span_literal") or ""
    target = next(
        (
            term.get("sub_span")
            for term in claim.get("terms", [])
            if term.get("term_type") == "celltype" and term.get("sub_span")
        ),
        None,
    )
    target_position = span.find(target) if target else None
    if target_position is not None and target_position < 0:
        target_position = None

    anchors: list[tuple[int, str]] = []
    seen: set[str] = set()
    for term in claim.get("terms", []):
        sub_span = term.get("sub_span")
        if not sub_span or sub_span in seen:
            continue
        seen.add(sub_span)
        anchors.append((_term_position(span, sub_span, target_position), sub_span))
    anchors.sort()
    if not anchors or not any(
        term.get("term_type") == "gene" for term in claim.get("terms", [])
    ):
        return None

    candidates: list[str] = []
    first = anchors[0][1]
    for start in _all_offsets(manuscript_text, first):
        cursor = start + len(first)
        end = cursor
        matched = True
        for _, anchor in anchors[1:]:
            position = manuscript_text.find(anchor, cursor, start + 2_000)
            if position < 0:
                matched = False
                break
            cursor = position + len(anchor)
            end = cursor
        if matched:
            candidates.append(manuscript_text[start:end])
    if not candidates:
        return None
    candidates.sort(key=len)
    if len(candidates) > 1 and len(candidates[0]) == len(candidates[1]):
        return None
    return candidates[0]


def _expand_gene_shorthand(span: str, label: str) -> str | None:
    """Return an exact shared source token such as FOXP1/2 for FOXP2."""

    match = re.fullmatch(r"(.+?)(\d+)", label)
    if match is None:
        return None
    prefix, suffix = match.groups()
    pattern = re.compile(
        rf"\b{re.escape(prefix)}(?:\d+/{re.escape(suffix)}|{re.escape(suffix)}/\d+)\b"
    )
    shorthand = pattern.search(span)
    return shorthand.group(0) if shorthand else None


def prepare_raw_claims(
    manuscript_text: str, raw_claims: Iterable[dict[str, Any]]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Align source spans and exclude records that are not marker claims."""

    prepared: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    reanchored = 0
    reconstructed = 0
    expanded_shorthand = 0
    implicit_unaligned_terms = 0
    excluded_terms: list[dict[str, Any]] = []
    raw_list = list(raw_claims)
    for index, raw in enumerate(raw_list):
        claim = copy.deepcopy(raw)
        genes = [
            term
            for term in claim.get("terms", [])
            if term.get("term_type") == "gene"
        ]
        if not genes:
            excluded.append({"raw_index": index, "reason": "no_marker_gene"})
            continue
        explicit_genes = [gene for gene in genes if gene.get("sub_span")]
        for gene in genes:
            if not gene.get("sub_span"):
                excluded_terms.append(
                    {
                        "raw_index": index,
                        "normalized_label": gene.get("normalized_label"),
                        "reason": "implicit_marker_gene",
                    }
                )
        if not explicit_genes:
            excluded.append(
                {"raw_index": index, "reason": "no_explicit_marker_gene"}
            )
            continue
        claim["terms"] = [
            term
            for term in claim.get("terms", [])
            if term.get("term_type") != "gene" or term.get("sub_span")
        ]

        span = claim.get("span_literal") or ""
        exact_span = _reanchor_span(manuscript_text, span)
        if exact_span is None:
            exact_span = _reconstruct_span_from_terms(manuscript_text, claim)
            if exact_span is not None:
                reconstructed += 1
        if exact_span is not None and exact_span != span:
            claim["span_literal"] = exact_span
            reanchored += 1
        aligned_span = claim.get("span_literal") or ""
        for term in claim.get("terms", []):
            sub_span = term.get("sub_span")
            if (
                term.get("term_type") == "gene"
                and sub_span
                and sub_span not in aligned_span
            ):
                shorthand = _expand_gene_shorthand(
                    aligned_span, term.get("normalized_label") or ""
                )
                if shorthand:
                    term["sub_span"] = shorthand
                    expanded_shorthand += 1
            elif sub_span and sub_span not in aligned_span:
                term["sub_span"] = None
                implicit_unaligned_terms += 1
        prepared.append(claim)

    report = {
        "raw_claims": len(raw_list),
        "retained_claims": len(prepared),
        "reanchored_spans": reanchored,
        "reconstructed_spans": reconstructed,
        "expanded_gene_shorthand": expanded_shorthand,
        "implicit_unaligned_terms": implicit_unaligned_terms,
        "excluded_claims": excluded,
        "excluded_terms": excluded_terms,
    }
    return prepared, report


def _claim_id(source_hash: str, span_literal: str, target_label: str) -> str:
    key = f"{source_hash}\t{span_literal}\t{target_label}"
    return "claim:" + hashlib.sha256(key.encode("utf-8")).hexdigest()[:20]


def _coalesce_claims(claims: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Merge records for the same source span and target into one marker panel."""

    output: list[dict[str, Any]] = []
    by_id: dict[str, dict[str, Any]] = {}
    for claim in claims:
        claim_id = claim["claim_id"]
        existing = by_id.get(claim_id)
        if existing is None:
            output.append(claim)
            by_id[claim_id] = claim
            continue
        if (
            existing["span_literal"] != claim["span_literal"]
            or existing["span_offset"] != claim["span_offset"]
        ):
            output.append(claim)
            continue
        if claim["summary"] not in existing["summary"]:
            existing["summary"] = f"{existing['summary']} {claim['summary']}"
        term_fingerprints = {json_fingerprint(term) for term in existing["terms"]}
        for term in claim["terms"]:
            fingerprint = json_fingerprint(term)
            if fingerprint not in term_fingerprints:
                existing["terms"].append(term)
                term_fingerprints.add(fingerprint)
    return output


def json_fingerprint(value: Any) -> str:
    """Return a stable digest for a JSON-compatible value."""

    import json

    payload = json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def make_claim_document(
    *, source_id: str, manuscript_text: str, raw_claims: Iterable[dict[str, Any]]
) -> dict[str, Any]:
    """Add deterministic provenance fields to raw LLM claim objects."""

    source_hash = sha256_text(manuscript_text)
    claims: list[dict[str, Any]] = []
    for raw in raw_claims:
        span = raw.get("span_literal") or ""
        span_offset = _find_offset(manuscript_text, span)
        terms: list[dict[str, Any]] = []
        for raw_term in raw.get("terms", []):
            sub_span = raw_term.get("sub_span")
            term = {
                "sub_span": sub_span,
                "sub_offset": _find_offset(span, sub_span),
                "normalized_label": raw_term.get("normalized_label"),
                "term_type": raw_term.get("term_type"),
                "provenance": "explicit" if sub_span else "implicit",
            }
            if raw_term.get("term_type") == "gene":
                term["direction"] = raw_term.get("direction")
            terms.append(term)

        target_labels = [
            term.get("normalized_label") or ""
            for term in terms
            if term.get("term_type") == "celltype"
        ]
        target_label = target_labels[0] if len(target_labels) == 1 else ""
        if manuscript_text.count(span) > 1 and any(
            term["term_type"] == "celltype" and term["provenance"] == "implicit"
            for term in terms
        ):
            span_offset = None
        claims.append(
            {
                "claim_id": _claim_id(source_hash, span, target_label),
                "span_literal": span,
                "span_offset": span_offset,
                "summary": raw.get("summary") or "",
                "terms": terms,
            }
        )

    return {
        "schema_version": CLAIMS_SCHEMA,
        "source": {"id": source_id, "sha256": source_hash},
        "producer": {"name": "mrkr", "version": __version__},
        "claims": _coalesce_claims(claims),
    }


def _offset_matches(container: str, value: str, offset: Any) -> bool:
    return (
        isinstance(offset, list)
        and len(offset) == 2
        and all(isinstance(item, int) for item in offset)
        and 0 <= offset[0] <= offset[1] <= len(container)
        and container[offset[0] : offset[1]] == value
    )


def validate_document(
    document: Any, manuscript_text: str | None = None
) -> dict[str, Any]:
    """Validate schema, biological cardinality, and exact source alignment."""

    errors: list[ValidationIssue] = []
    warnings: list[ValidationIssue] = []

    def add(code: str, path: str, message: str) -> None:
        errors.append(ValidationIssue(code, path, message))

    def warn(code: str, path: str, message: str) -> None:
        warnings.append(ValidationIssue(code, path, message))

    if not isinstance(document, dict):
        add("document.type", "$", "document must be a JSON object")
        return _report(0, errors, warnings)

    schema = document.get("schema_version")
    if schema not in {CLAIMS_SCHEMA, ONTO_SCHEMA}:
        add("schema.unsupported", "schema_version", f"expected {CLAIMS_SCHEMA} or {ONTO_SCHEMA}")
    ols_results: set[tuple[str, str | None, bool | None]] = set()
    if schema == CLAIMS_SCHEMA and "grounding" in document:
        add(
            "grounding.premature",
            "grounding",
            "extraction documents cannot contain grounding metadata",
        )
    if schema == ONTO_SCHEMA:
        grounding = document.get("grounding")
        if not isinstance(grounding, dict):
            add("grounding.missing", "grounding", "onto documents require grounding metadata")
            grounding = {}
        genes = grounding.get("genes")
        if not isinstance(genes, dict):
            add("grounding.genes", "grounding.genes", "gene-map metadata is required")
        else:
            if genes.get("provider") != "offline-gene-map":
                add("grounding.genes", "grounding.genes.provider", "unsupported gene-map provider")
            if not isinstance(genes.get("organism"), str) or not genes.get("organism"):
                add("grounding.organism", "grounding.genes.organism", "organism is required")
            digest = genes.get("sha256")
            if not isinstance(digest, str) or not digest.startswith("sha256:"):
                add("grounding.genes", "grounding.genes.sha256", "gene-map SHA-256 is required")
            canonical_digest = genes.get("canonical_sha256")
            if canonical_digest is not None and (
                not isinstance(canonical_digest, str)
                or not canonical_digest.startswith("sha256:")
            ):
                add(
                    "grounding.genes",
                    "grounding.genes.canonical_sha256",
                    "canonical gene-map digest must be a SHA-256",
                )
        organism_grounding = grounding.get("organism")
        if not isinstance(organism_grounding, dict):
            add(
                "grounding.organism",
                "grounding.organism",
                "organism grounding metadata is required",
            )
            organism_grounding = {}
        if organism_grounding.get("provider") != "NCBI Taxonomy":
            add(
                "grounding.organism",
                "grounding.organism.provider",
                "organism grounding must use NCBI Taxonomy",
            )
        if not isinstance(organism_grounding.get("label"), str) or not organism_grounding.get(
            "label"
        ):
            add(
                "grounding.organism",
                "grounding.organism.label",
                "organism label is required",
            )
        organism_identifier = organism_grounding.get("ontology_term")
        if not isinstance(organism_identifier, str) or not organism_identifier.startswith(
            "NCBITaxon:"
        ):
            add(
                "grounding.organism",
                "grounding.organism.ontology_term",
                "NCBI Taxonomy identifier is required",
            )
        service = grounding.get("ontology_service")
        if not isinstance(service, dict):
            add("grounding.service", "grounding.ontology_service", "ontology-service metadata is required")
            queries: list[Any] = []
        else:
            if service.get("provider") != "OLS4":
                add("grounding.service", "grounding.ontology_service.provider", "unsupported ontology service")
            if not isinstance(service.get("endpoint"), str) or not service.get("endpoint"):
                add("grounding.service", "grounding.ontology_service.endpoint", "service endpoint is required")
            queries = service.get("queries")
            if not isinstance(queries, list):
                add("grounding.queries", "grounding.ontology_service.queries", "queries must be an array")
                queries = []
        for query_index, query in enumerate(queries):
            query_path = f"grounding.ontology_service.queries[{query_index}]"
            if not isinstance(query, dict):
                add("grounding.query", query_path, "query evidence must be an object")
                continue
            ontology = query.get("ontology")
            ontology_term = query.get("ontology_term")
            exact = query.get("exact")
            if not isinstance(query.get("query"), str) or not query.get("query"):
                add("grounding.query", f"{query_path}.query", "query text is required")
            if ontology not in {"cl", "uberon"}:
                add("grounding.query", f"{query_path}.ontology", "ontology must be cl or uberon")
            if not isinstance(query.get("retrieved_at"), str) or not query.get("retrieved_at"):
                add("grounding.query", f"{query_path}.retrieved_at", "retrieval time is required")
            response_digest = query.get("response_sha256")
            if not isinstance(response_digest, str) or not response_digest.startswith("sha256:"):
                add("grounding.query", f"{query_path}.response_sha256", "response SHA-256 is required")
            if ontology_term is None and exact is not None:
                add("grounding.query", query_path, "unresolved query must have exact=null")
            elif ontology_term is not None and not isinstance(exact, bool):
                add("grounding.query", query_path, "resolved query requires Boolean exactness")
            if ontology in {"cl", "uberon"}:
                ols_results.add((ontology, ontology_term, exact))

    source = document.get("source")
    if not isinstance(source, dict):
        add("source.missing", "source", "source metadata is required")
        source = {}
    if not source.get("id"):
        add("source.id", "source.id", "source id is required")
    if not source.get("sha256"):
        add("source.sha256", "source.sha256", "source SHA-256 digest is required")
    if manuscript_text is not None and source.get("sha256") != sha256_text(manuscript_text):
        add("source.hash_mismatch", "source.sha256", "digest does not match the manuscript")

    claims = document.get("claims")
    if not isinstance(claims, list):
        add("claims.type", "claims", "claims must be a JSON array")
        return _report(0, errors, warnings)

    seen_claim_ids: set[str] = set()
    for claim_index, claim in enumerate(claims):
        base = f"claims[{claim_index}]"
        if not isinstance(claim, dict):
            add("claim.type", base, "claim must be a JSON object")
            continue

        claim_id = claim.get("claim_id")
        if not isinstance(claim_id, str) or not claim_id:
            add("claim.id", f"{base}.claim_id", "claim id is required")
        elif claim_id in seen_claim_ids:
            add("claim.id_duplicate", f"{base}.claim_id", "claim id must be unique")
        else:
            seen_claim_ids.add(claim_id)

        span = claim.get("span_literal")
        summary = claim.get("summary")
        if not isinstance(span, str) or not span:
            add("claim.span", f"{base}.span_literal", "verbatim evidence span is required")
            span = ""
        if not isinstance(summary, str) or not summary:
            add("claim.summary", f"{base}.summary", "normalized summary is required")
            summary = ""

        if manuscript_text is not None and not _offset_matches(
            manuscript_text, span, claim.get("span_offset")
        ):
            add(
                "claim.span_offset",
                f"{base}.span_offset",
                "offset must select span_literal exactly from the manuscript",
            )

        terms = claim.get("terms")
        if not isinstance(terms, list):
            add("terms.type", f"{base}.terms", "terms must be a JSON array")
            continue

        celltype_count = 0
        gene_count = 0
        organism_count = 0
        for term_index, term in enumerate(terms):
            term_path = f"{base}.terms[{term_index}]"
            if not isinstance(term, dict):
                add("term.type", term_path, "term must be a JSON object")
                continue
            term_type = term.get("term_type")
            if term_type not in TERM_TYPES:
                add("term.term_type", f"{term_path}.term_type", "unsupported term type")
            celltype_count += term_type == "celltype"
            gene_count += term_type == "gene"
            organism_count += term_type == "organism"

            label = term.get("normalized_label")
            if not isinstance(label, str) or not label:
                add("term.label", f"{term_path}.normalized_label", "normalized label is required")
            elif isinstance(summary, str) and label.casefold() not in summary.casefold():
                report = add if term_type == "organism" else warn
                report(
                    "term.label_missing_summary",
                    f"{term_path}.normalized_label",
                    "normalized label should occur in the claim summary",
                )

            sub_span = term.get("sub_span")
            provenance = term.get("provenance")
            if sub_span is None:
                if term_type == "gene":
                    add(
                        "term.gene_explicit",
                        f"{term_path}.sub_span",
                        "marker genes must be explicit in the evidence span",
                    )
                if provenance != "implicit":
                    add("term.provenance", f"{term_path}.provenance", "missing sub_span requires implicit provenance")
                if term.get("sub_offset") is not None:
                    add("term.sub_offset", f"{term_path}.sub_offset", "implicit term cannot have a source offset")
            elif not isinstance(sub_span, str) or not sub_span:
                add("term.sub_span", f"{term_path}.sub_span", "sub_span must be non-empty or null")
            else:
                if provenance != "explicit":
                    add("term.provenance", f"{term_path}.provenance", "sub_span requires explicit provenance")
                if not _offset_matches(span, sub_span, term.get("sub_offset")):
                    add(
                        "term.sub_offset",
                        f"{term_path}.sub_offset",
                        "offset must select sub_span exactly from span_literal",
                    )

            if term_type == "gene":
                if term.get("direction") not in GENE_DIRECTIONS:
                    add("term.direction", f"{term_path}.direction", "gene direction must be positive or negative")
            elif term.get("direction") is not None:
                add("term.direction", f"{term_path}.direction", "direction is only valid for gene terms")

            has_ontology_term = "ontology_term" in term
            has_exact = "exact" in term
            ontology_term = term.get("ontology_term")
            exact = term.get("exact")
            if schema == CLAIMS_SCHEMA and (has_ontology_term or has_exact):
                add(
                    "term.premature_grounding",
                    term_path,
                    "extraction output cannot contain ontology fields",
                )
            if schema == ONTO_SCHEMA:
                if not has_ontology_term or not has_exact:
                    add(
                        "term.grounding_state",
                        term_path,
                        "onto terms require ontology_term and exact fields",
                    )
                elif ontology_term is None and exact is not None:
                    add("term.grounding_state", term_path, "unresolved term must have exact=null")
                elif ontology_term is not None and not isinstance(exact, bool):
                    add("term.grounding_state", term_path, "grounded term must have exact=true or exact=false")
                if term_type == "gene" and ontology_term is not None and not str(
                    ontology_term
                ).startswith("ENS"):
                    add(
                        "term.gene_identifier",
                        f"{term_path}.ontology_term",
                        "gene identifiers must be Ensembl identifiers",
                    )
                if term_type == "organism" and ontology_term is not None and not str(
                    ontology_term
                ).startswith("NCBITaxon:"):
                    add(
                        "term.organism_identifier",
                        f"{term_path}.ontology_term",
                        "organism identifiers must be NCBI Taxonomy identifiers",
                    )
                if term_type == "organism" and (
                    ontology_term != organism_grounding.get("ontology_term")
                    or label.casefold()
                    != str(organism_grounding.get("label", "")).casefold()
                ):
                    add(
                        "term.organism_mismatch",
                        term_path,
                        "organism term must match document grounding metadata",
                    )
                ontology = {
                    "celltype": "cl",
                    "comparison": "cl",
                    "tissue": "uberon",
                }.get(term_type)
                if ontology and (ontology, ontology_term, exact) not in ols_results:
                    add(
                        "term.grounding_evidence",
                        term_path,
                        "term grounding must match recorded ontology query evidence",
                    )
                if ontology_term is None:
                    warn(
                        "term.grounding_unresolved",
                        term_path,
                        "term is present but has no stable identifier",
                    )
                elif exact is False:
                    warn(
                        "term.grounding_coarse",
                        term_path,
                        "term is grounded to a broader concept",
                    )

        if celltype_count != 1:
            add("claim.celltype_count", f"{base}.terms", "claim must contain exactly one target cell type")
        if gene_count < 1:
            add("claim.gene_count", f"{base}.terms", "claim must contain at least one marker gene")
        if organism_count != 1:
            add(
                "claim.organism_count",
                f"{base}.terms",
                "claim must contain exactly one organism",
            )

    return _report(len(claims), errors, warnings)


def _report(
    n_claims: int,
    errors: list[ValidationIssue],
    warnings: list[ValidationIssue],
) -> dict[str, Any]:
    return {
        "valid": not errors,
        "n_claims": n_claims,
        "n_errors": len(errors),
        "errors": [issue.as_dict() for issue in errors],
        "n_warnings": len(warnings),
        "warnings": [issue.as_dict() for issue in warnings],
    }


def assert_valid_document(
    document: Any, manuscript_text: str | None = None
) -> dict[str, Any]:
    """Return a validation report or raise ClaimValidationError."""

    report = validate_document(document, manuscript_text)
    if not report["valid"]:
        raise ClaimValidationError(report)
    return report
