"""Grounded marker-profile query utilities."""

from __future__ import annotations

import json
import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import anthropic
from pydantic import BaseModel, Field

from .config import config
from .llm import extract_json_from_response
from .map import load_gene_map, resolve_gene_id


GENERIC_LABEL_TOKENS = {
    "A",
    "AN",
    "AND",
    "ARE",
    "AS",
    "ASKING",
    "CELL",
    "CELLS",
    "COMPARE",
    "CONTEXT",
    "DO",
    "DOES",
    "EXPERIMENTAL",
    "FOR",
    "FROM",
    "GENE",
    "GENES",
    "IN",
    "IS",
    "LIST",
    "MARKER",
    "MARKERS",
    "OF",
    "OR",
    "PAPER",
    "PAPERS",
    "PROFILE",
    "PROFILES",
    "QUERY",
    "SAME",
    "STUDIES",
    "STUDY",
    "THE",
    "THESE",
    "THING",
    "TO",
    "TYPE",
    "TYPES",
    "USING",
    "WHAT",
    "WHERE",
    "WHICH",
    "WITH",
}

LABEL_STOPWORDS = {
    "CELL",
    "CELLS",
    "TYPE",
    "TYPES",
    "POPULATION",
    "POPULATIONS",
    "SUBPOPULATION",
    "SUBPOPULATIONS",
    "REPORTED",
    "OTHER",
    "LIKE",
}


@dataclass(frozen=True)
class QueryMarker:
    label: str
    feature_id: str | None


@dataclass(frozen=True)
class StructuredQuery:
    raw_query: str
    cell_type_label: str
    context: str
    markers: tuple[QueryMarker, ...]
    parser: str = "heuristic"


@dataclass(frozen=True)
class Profile:
    profile_id: int
    paper_id: int
    collection: str
    organism: str
    group_name: str
    text_blob: str
    paper_context_blob: str
    gene_names: tuple[str, ...]
    gene_ids: tuple[str, ...]
    gene_labels_by_id: dict[str, tuple[str, ...]]
    evidence_sentences: tuple[str, ...]
    doi: str | None
    title: str | None
    year: int | None


class LLMQueryMarker(BaseModel):
    label: str = Field(description="Reported gene symbol or marker label from the user query")
    feature_id: str | None = Field(default=None, description="Ensembl gene ID if known")


class LLMStructuredQuery(BaseModel):
    cell_type_label: str = Field(default="", description="Cell type label or common cell population name")
    context: str = Field(default="", description="Biological or experimental context")
    markers: list[LLMQueryMarker] = Field(default_factory=list, description="Marker genes from the query")


def normalize_label(text: str) -> str:
    return " ".join(str(text or "").upper().split())


def tokenize(text: str) -> tuple[str, ...]:
    tokens = []
    for token in re.findall(r"[A-Za-z0-9]+", normalize_label(text)):
        if len(token) < 2 and token not in {"B", "T"}:
            continue
        if token in GENERIC_LABEL_TOKENS:
            continue
        tokens.append(token)
    return tuple(tokens)


def label_tokens(text: str) -> set[str]:
    return {token for token in tokenize(text) if token not in LABEL_STOPWORDS}


def token_jaccard(left: Iterable[str], right: Iterable[str]) -> float:
    left_set = set(left)
    right_set = set(right)
    union = left_set | right_set
    if not union:
        return 0.0
    return len(left_set & right_set) / len(union)


def label_similarity(query_label: str, profile_label: str) -> tuple[str, float]:
    query_norm = normalize_label(query_label)
    profile_norm = normalize_label(profile_label)
    if not query_norm:
        return "none", 0.0
    if query_norm == profile_norm:
        return "exact", 1.0
    if query_norm in profile_norm or profile_norm in query_norm:
        return "partial", 0.75

    score = token_jaccard(label_tokens(query_norm), label_tokens(profile_norm))
    if score > 0:
        return "partial", score
    return "different", 0.0


def extract_candidate_gene_labels(text: str, gene_map: dict[str, str]) -> list[QueryMarker]:
    candidates: list[QueryMarker] = []
    seen: set[str] = set()
    for token in re.findall(r"[A-Za-z0-9][A-Za-z0-9\-]*", text):
        cleaned = token.rstrip("+").strip()
        if len(cleaned) < 3 and cleaned.upper() not in {"CD4", "CD8"}:
            continue
        feature_id = resolve_gene_id(cleaned, gene_map)
        if not feature_id:
            continue
        key = feature_id
        if key in seen:
            continue
        seen.add(key)
        candidates.append(QueryMarker(label=cleaned.upper(), feature_id=feature_id))
    return candidates


def infer_context(raw_query: str) -> str:
    match = re.search(r"\b(?:in|from|within|among)\s+(.+)$", raw_query, flags=re.IGNORECASE)
    if not match:
        return ""
    context = re.sub(r"[?.!]\s*$", "", match.group(1)).strip()
    return context


def infer_label(raw_query: str, markers: tuple[QueryMarker, ...], context: str) -> str:
    text = raw_query
    if context:
        text = re.sub(r"\b(?:in|from|within|among)\s+" + re.escape(context) + r"\s*$", "", text, flags=re.IGNORECASE)
    for marker in markers:
        text = re.sub(r"\b" + re.escape(marker.label) + r"\+?\b", " ", text, flags=re.IGNORECASE)
    tokens = [token for token in tokenize(text) if token not in {marker.label for marker in markers}]
    if not tokens:
        return ""

    # Keep the biologically meaningful tail of verbose questions.
    if len(tokens) > 6:
        tokens = tokens[-6:]
    return " ".join(tokens)


def parse_query_heuristic(
    raw_query: str,
    *,
    label: str | None = None,
    context: str | None = None,
    gene_labels: tuple[str, ...] = (),
    gene_map_file: Path | None = None,
) -> StructuredQuery:
    gene_map = load_gene_map(gene_map_file=gene_map_file)
    marker_by_id: dict[str, QueryMarker] = {}

    for marker in extract_candidate_gene_labels(raw_query, gene_map):
        if marker.feature_id:
            marker_by_id[marker.feature_id] = marker

    for gene_label in gene_labels:
        feature_id = resolve_gene_id(gene_label, gene_map)
        key = feature_id or normalize_label(gene_label)
        marker_by_id[key] = QueryMarker(label=normalize_label(gene_label), feature_id=feature_id)

    markers = tuple(marker_by_id.values())
    parsed_context = context if context is not None else infer_context(raw_query)
    parsed_label = label if label is not None else infer_label(raw_query, markers, parsed_context)

    return StructuredQuery(
        raw_query=raw_query,
        cell_type_label=parsed_label.strip(),
        context=parsed_context.strip(),
        markers=markers,
        parser="heuristic",
    )


def parse_query_llm(
    raw_query: str,
    *,
    label: str | None = None,
    context: str | None = None,
    gene_labels: tuple[str, ...] = (),
    gene_map_file: Path | None = None,
) -> StructuredQuery:
    config.validate()
    client = anthropic.Anthropic(api_key=config.anthropic_api_key)
    prompt = f"""
Parse this biological marker-profile query into JSON for deterministic search.

Return only this JSON object:
{{
  "cell_type_label": "cell type or population label, empty string if absent",
  "context": "biological, disease, tissue, assay, or experimental context, empty string if absent",
  "markers": [
    {{"label": "gene symbol or marker label", "feature_id": null}}
  ]
}}

Rules:
- Do not answer the biological question.
- Include only markers explicitly named by the user or supplied in gene overrides.
- Preserve marker labels as gene symbols or surface markers where they appear.
- A pattern such as "CCR8+ Tregs" means cell_type_label="Treg" and markers=[{{"label":"CCR8","feature_id":null}}].
- A phrase such as "TREM2+ macrophages in tumors" means cell_type_label="macrophage", context="tumors", and markers=[{{"label":"TREM2","feature_id":null}}].
- Put tissue, disease, cancer type, technology, and experimental setting in context.
- If the user provides overrides, respect them.
- Prefer concise labels such as "Treg", "macrophage", "exhausted T cell", "monocyte", or "CD8 T cell".

User query: {raw_query}
Label override: {label or ""}
Context override: {context or ""}
Gene overrides: {", ".join(gene_labels)}
"""
    message = client.messages.create(
        model=config.anthropic_model,
        max_tokens=1000,
        temperature=0.0,
        messages=[{"role": "user", "content": prompt}],
    )
    parsed = LLMStructuredQuery.model_validate_json(extract_json_from_response(message.content[0].text))
    gene_map = load_gene_map(gene_map_file=gene_map_file)
    markers = []
    seen = set()
    for marker in parsed.markers:
        feature_id = marker.feature_id or resolve_gene_id(marker.label, gene_map)
        key = feature_id or normalize_label(marker.label)
        if key in seen:
            continue
        seen.add(key)
        markers.append(QueryMarker(label=normalize_label(marker.label), feature_id=feature_id))
    for gene_label in gene_labels:
        feature_id = resolve_gene_id(gene_label, gene_map)
        key = feature_id or normalize_label(gene_label)
        if key in seen:
            continue
        seen.add(key)
        markers.append(QueryMarker(label=normalize_label(gene_label), feature_id=feature_id))
    return StructuredQuery(
        raw_query=raw_query,
        cell_type_label=(label if label is not None else parsed.cell_type_label).strip(),
        context=(context if context is not None else parsed.context).strip(),
        markers=tuple(markers),
        parser="llm",
    )


def load_profiles(db_path: Path, organism: str | None = "homo_sapiens") -> list[Profile]:
    if not db_path.exists():
        raise FileNotFoundError(f"LLMarkers SQLite database not found: {db_path}")
    query = """
        SELECT
            p.profile_id,
            p.paper_id,
            p.collection,
            p.organism,
            p.group_name,
            p.text_blob,
            p.paper_context_blob,
            p.gene_names_json,
            p.gene_ids_json,
            p.evidence_sentences_json,
            pa.doi,
            pa.title,
            pa.year
        FROM profiles AS p
        JOIN papers AS pa ON pa.paper_id = p.paper_id
        WHERE p.n_gene_ids > 0
    """
    params: list[str] = []
    if organism:
        query += " AND p.organism = ?"
        params.append(organism)
    query += " ORDER BY p.profile_id"

    with sqlite3.connect(db_path) as conn:
        rows = conn.execute(query, params).fetchall()
        marker_rows = conn.execute(
            """
            SELECT
                p.profile_id,
                m.feature_id,
                m.feature_name
            FROM profiles AS p
            JOIN markers AS m
              ON m.paper_id = p.paper_id
             AND COALESCE(m.organism, '') = COALESCE(p.organism, '')
             AND m.group_name = p.group_name
            WHERE COALESCE(m.feature_id, '') <> ''
            """
        ).fetchall()

    labels_by_profile_and_id: dict[int, dict[str, set[str]]] = {}
    for profile_id, feature_id, feature_name in marker_rows:
        labels_by_profile_and_id.setdefault(int(profile_id), {}).setdefault(str(feature_id), set()).add(str(feature_name))

    profiles = []
    for row in rows:
        profiles.append(
            Profile(
                profile_id=int(row[0]),
                paper_id=int(row[1]),
                collection=str(row[2] or ""),
                organism=str(row[3] or ""),
                group_name=str(row[4] or ""),
                text_blob=str(row[5] or ""),
                paper_context_blob=str(row[6] or ""),
                gene_names=tuple(json.loads(row[7] or "[]")),
                gene_ids=tuple(json.loads(row[8] or "[]")),
                gene_labels_by_id={
                    gene_id: tuple(sorted(labels))
                    for gene_id, labels in labels_by_profile_and_id.get(int(row[0]), {}).items()
                },
                evidence_sentences=tuple(json.loads(row[9] or "[]")),
                doi=row[10],
                title=row[11],
                year=row[12],
            )
        )
    return profiles


def marker_similarity(query_gene_ids: set[str], profile: Profile) -> dict[str, object]:
    profile_gene_ids = set(profile.gene_ids)
    shared_ids = query_gene_ids & profile_gene_ids
    union = query_gene_ids | profile_gene_ids

    def labels_for(gene_id: str) -> str:
        labels = profile.gene_labels_by_id.get(gene_id)
        if labels:
            return "/".join(labels)
        return gene_id

    shared_genes = [labels_for(gene_id) for gene_id in sorted(shared_ids)]
    profile_only_genes = [labels_for(gene_id) for gene_id in sorted(profile_gene_ids - query_gene_ids)]
    return {
        "shared_gene_ids": sorted(shared_ids),
        "shared_genes": shared_genes,
        "profile_only_genes": profile_only_genes,
        "query_coverage": len(shared_ids) / len(query_gene_ids) if query_gene_ids else None,
        "profile_coverage": len(shared_ids) / len(profile_gene_ids) if profile_gene_ids else 0.0,
        "jaccard": len(shared_ids) / len(union) if union else None,
    }


def rank_positions(scored: list[tuple[int, float]], *, limit: int = 200) -> dict[int, int]:
    ranked = sorted(scored, key=lambda item: (-item[1], item[0]))
    return {profile_id: rank for rank, (profile_id, score) in enumerate(ranked[:limit], start=1) if score > 0}


def evidence_sentence(profile: Profile, shared_genes: list[str]) -> str:
    if not profile.evidence_sentences:
        return ""
    if shared_genes:
        for sentence in profile.evidence_sentences:
            sentence_norm = normalize_label(sentence)
            if any(normalize_label(gene) in sentence_norm for gene in shared_genes):
                return sentence
    return profile.evidence_sentences[0]


def interpret_match(
    query: StructuredQuery,
    marker: dict[str, object],
    label_relation_name: str,
    context_score: float,
) -> str:
    query_coverage = marker["query_coverage"]
    jaccard = marker["jaccard"]
    if query_coverage is None:
        if label_relation_name in {"exact", "partial"} and context_score > 0:
            return "label and context match"
        if label_relation_name in {"exact", "partial"}:
            return "label match"
        return "context/text match"

    if query_coverage == 1.0 and (jaccard or 0) >= 0.5 and label_relation_name in {"exact", "partial"}:
        return "strong marker-program match with similar label"
    if query_coverage == 1.0 and (jaccard or 0) >= 0.5:
        return "possible alias or cross-label marker-program match"
    if query_coverage == 1.0:
        return "query marker pattern present in broader profile"
    if query_coverage > 0 and label_relation_name in {"exact", "partial"}:
        return "similar label with partial marker overlap"
    if query_coverage > 0:
        return "partial marker-program match"
    if label_relation_name in {"exact", "partial"}:
        return "similar label without queried marker overlap"
    return "weak deterministic match"


def resolve_query(
    db_path: Path,
    query: StructuredQuery,
    *,
    organism: str | None = "homo_sapiens",
    top_n: int = 10,
) -> dict[str, object]:
    profiles = load_profiles(db_path, organism=organism)
    query_gene_ids = {marker.feature_id for marker in query.markers if marker.feature_id}
    query_context_tokens = tokenize(query.context)
    query_text_tokens = tokenize(" ".join(part for part in [query.cell_type_label, query.context] if part))

    profile_scores: dict[int, dict[str, object]] = {}
    marker_scored: list[tuple[int, float]] = []
    label_scored: list[tuple[int, float]] = []
    context_scored: list[tuple[int, float]] = []

    for profile in profiles:
        marker = marker_similarity(query_gene_ids, profile)
        label_relation_name, label_score = label_similarity(query.cell_type_label, profile.group_name)
        claim_text_score = token_jaccard(query_text_tokens, tokenize(profile.text_blob))
        context_text_score = token_jaccard(query_context_tokens, tokenize(profile.paper_context_blob)) if query_context_tokens else 0.0
        context_score = max(claim_text_score, context_text_score)

        marker_rank_score = 0.0
        if query_gene_ids:
            marker_rank_score = (
                2.0 * float(marker["query_coverage"] or 0)
                + float(marker["jaccard"] or 0)
                + 0.25 * float(marker["profile_coverage"] or 0)
            )

        profile_scores[profile.profile_id] = {
            "profile": profile,
            "marker": marker,
            "label_relation": label_relation_name,
            "label_score": label_score,
            "claim_text_score": claim_text_score,
            "context_score": context_score,
            "marker_rank_score": marker_rank_score,
        }
        marker_scored.append((profile.profile_id, marker_rank_score))
        label_scored.append((profile.profile_id, label_score))
        context_scored.append((profile.profile_id, context_score))

    ranks = {
        "marker_nearest": rank_positions(marker_scored),
        "label_nearest": rank_positions(label_scored),
        "context_nearest": rank_positions(context_scored),
    }
    consensus_scores: list[tuple[int, float]] = []
    for profile_id in profile_scores:
        score = 0.0
        for rank_map in ranks.values():
            rank = rank_map.get(profile_id)
            if rank is not None:
                score += 1.0 / (60 + rank)
        if score > 0:
            consensus_scores.append((profile_id, score))
    ranks["consensus"] = rank_positions(consensus_scores)

    candidate_ids = set()
    if query_gene_ids:
        # If the user supplied markers, matches with no shared queried genes are
        # usually not useful as primary candidates. Label/context-only profiles
        # can still be inspected through the ranked_views block in JSON output.
        candidate_ids.update(
            profile_id
            for profile_id, rank in ranks["marker_nearest"].items()
            if rank <= max(top_n * 5, 25)
        )
    else:
        for rank_map in ranks.values():
            candidate_ids.update(profile_id for profile_id, rank in rank_map.items() if rank <= max(top_n, 25))

    used_label_context_fallback = False
    if not candidate_ids:
        used_label_context_fallback = True
        for rank_map in [ranks["label_nearest"], ranks["context_nearest"]]:
            candidate_ids.update(profile_id for profile_id, rank in rank_map.items() if rank <= max(top_n, 25))

    consensus_by_profile = dict(consensus_scores)

    def match_sort_key(profile_id: int) -> tuple[float, ...]:
        row = profile_scores[profile_id]
        marker = row["marker"]
        if query_gene_ids:
            return (
                -float(marker["query_coverage"] or 0.0),
                -float(row["label_score"] or 0.0),
                -float(marker["jaccard"] or 0.0),
                -float(marker["profile_coverage"] or 0.0),
                -float(row["context_score"] or 0.0),
                -consensus_by_profile.get(profile_id, 0.0),
                float(profile_id),
            )
        return (
            -float(row["label_score"] or 0.0),
            -float(row["context_score"] or 0.0),
            -consensus_by_profile.get(profile_id, 0.0),
            float(profile_id),
        )

    matches = []
    for profile_id in sorted(candidate_ids, key=match_sort_key)[:top_n]:
        row = profile_scores[profile_id]
        profile = row["profile"]
        marker = row["marker"]
        shared_genes = marker["shared_genes"]
        match = {
            "profile_id": profile.profile_id,
            "label": profile.group_name,
            "paper": {
                "paper_id": profile.paper_id,
                "title": profile.title,
                "doi": profile.doi,
                "year": profile.year,
                "collection": profile.collection,
            },
            "ranks": {
                view: rank_map.get(profile.profile_id)
                for view, rank_map in ranks.items()
            },
            "marker_similarity": marker,
            "label_relation": row["label_relation"],
            "label_score": row["label_score"],
            "claim_text_score": row["claim_text_score"],
            "context_score": row["context_score"],
            "source_sentence": evidence_sentence(profile, shared_genes),
            "profile_genes": [
                "/".join(profile.gene_labels_by_id.get(gene_id, (gene_id,)))
                for gene_id in sorted(set(profile.gene_ids))
            ],
            "interpretation": interpret_match(
                query,
                marker,
                str(row["label_relation"]),
                float(row["context_score"]),
            ),
        }
        matches.append(match)

    return {
        "query": {
            "raw_query": query.raw_query,
            "parser": query.parser,
            "cell_type_label": query.cell_type_label,
            "context": query.context,
            "markers": [
                {"label": marker.label, "feature_id": marker.feature_id}
                for marker in query.markers
            ],
        },
        "database": str(db_path),
        "organism": organism,
        "n_profiles_searched": len(profiles),
        "notes": (
            [
                "No searched profile contained the resolved query marker genes; returned label/context matches instead."
            ]
            if used_label_context_fallback and query_gene_ids
            else []
        ),
        "ranked_views": {
            view: [
                {"profile_id": profile_id, "rank": rank}
                for profile_id, rank in sorted(rank_map.items(), key=lambda item: item[1])[:top_n]
            ]
            for view, rank_map in ranks.items()
        },
        "matches": matches,
    }


def print_query_result(result: dict[str, object], top_n: int) -> None:
    query = result["query"]
    print("\nStructured query")
    print("----------------")
    print(f"parser:  {query.get('parser') or '(unknown)'}")
    print(f"label:   {query['cell_type_label'] or '(none)'}")
    print(f"context: {query['context'] or '(none)'}")
    markers = query["markers"]
    if markers:
        marker_text = ", ".join(
            f"{marker['label']} ({marker['feature_id'] or 'unmapped'})"
            for marker in markers
        )
    else:
        marker_text = "(none)"
    print(f"markers: {marker_text}")
    print(f"\nSearched {result['n_profiles_searched']:,} profiles\n")
    for note in result.get("notes", []):
        print(f"Note: {note}\n")

    print(f"Top {top_n} consensus matches")
    print("------------------------")
    for i, match in enumerate(result["matches"], start=1):
        marker = match["marker_similarity"]
        paper = match["paper"]
        shared = ", ".join(marker["shared_genes"][:8]) or "-"
        profile_only = ", ".join(marker["profile_only_genes"][:8]) or "-"
        jaccard = marker["jaccard"]
        query_coverage = marker["query_coverage"]
        ranks = match["ranks"]
        print(f"{i}. profile {match['profile_id']} | {match['label']} | {paper['year'] or ''}")
        print(f"   paper: {paper['title'] or ''}")
        print(f"   doi:   {paper['doi'] or ''}")
        jaccard_text = f"{jaccard:.3f}" if isinstance(jaccard, float) else "NA"
        query_coverage_text = f"{query_coverage:.3f}" if isinstance(query_coverage, float) else "NA"
        print(
            f"   marker: J={jaccard_text} query_cov={query_coverage_text} "
            f"shared=[{shared}] profile_only=[{profile_only}]"
        )
        print(
            "   relation: "
            f"{match['label_relation']} | "
            f"label={match['label_score']:.3f} "
            f"context={match['context_score']:.3f} | "
            f"ranks marker/label/context/consensus="
            f"{ranks['marker_nearest']}/{ranks['label_nearest']}/{ranks['context_nearest']}/{ranks['consensus']}"
        )
        print(f"   call: {match['interpretation']}")
        if match["source_sentence"]:
            print(f"   evidence: {match['source_sentence'][:260]}")
        print()
