"""Programmatic grounding for mrkr claim objects.

The LLM (extract) emits spans + normalized labels + term_type (+ direction on genes) and NEVER an
ontology id. This module assigns the ids:

  gene                 -> Ensembl gene id      (offline gmap.txt, via mrkr.map)
  celltype, comparison -> Cell Ontology (CL)   (OLS tag_text)
  tissue               -> UBERON               (OLS tag_text)

CL/UBERON grounding rule (one span-coverage rule): singularize the query (CL stores singulars),
run OLS tag_text, take the term covering the LONGEST span of the label; full-label coverage ->
exact=True, sub-phrase (a base concept) -> exact=False (coarse); the universal root CL:0000000 is
blocklisted; nothing -> None. The stored normalized_label is never changed (query-time transform).

Genes stay on the offline gmap.txt map for reproducibility (no network); only CL/UBERON hit OLS.
"""

import copy
import hashlib
import json
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from .claims import ONTO_SCHEMA, assert_valid_document
from .map import load_gene_map, resolve_gene_id

_DELIM = " ,.;:!?\t\n()[]{}\"'/\\-"
_OLS_URL = "https://www.ebi.ac.uk/ols4/api/v2/tag_text"
_PACKAGED_GENE_MAP_ORGANISM = "homo_sapiens"
_ols_cache: dict = {}
_ols_evidence: dict = {}


class GroundingServiceError(RuntimeError):
    """Raised when ontology grounding could not be completed reliably."""


class GeneGroundingConflict(ValueError):
    """Raised when source and normalized gene labels resolve differently."""


def _iri2curie(iri: str) -> str:
    frag = iri.rstrip("/").split("/")[-1]
    return frag.replace("_", ":", 1) if "_" in frag else frag


def _singularize(label: str) -> str:
    def sw(w: str) -> str:
        lw = w.lower()
        if len(w) > 3 and lw.endswith("ies"):
            return w[:-3] + "y"
        if len(w) > 3 and lw.endswith("s") and not lw.endswith("ss"):
            return w[:-1]
        return w

    return " ".join(sw(w) for w in label.split())


def _ols_once(label: str, ont: str):
    key = (label.lower(), ont)
    if key in _ols_cache:
        return _ols_cache[key]
    params = [("minLength", "4"), ("includeSubstrings", "false"),
              ("includeObsoleteEntities", "false"), ("delimiters", _DELIM), ("ontologyId", ont)]
    url = _OLS_URL + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, data=json.dumps({"text": label}).encode(),
                                 headers={"Content-Type": "application/json"})
    res = (None, None)
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            d = json.load(r)
    except Exception as exc:
        raise GroundingServiceError(
            f"OLS request failed for {ont}:{label!r}: {exc}"
        ) from exc
    cand = []
    for e in (d.get("entities") or []):
        curie = _iri2curie(e.get("term_iri") or e.get("iri") or "")
        if curie == "CL:0000000":            # blocklist the universal root
            continue
        s, en = e.get("start", 0), e.get("end", 0)
        cand.append((en - s, s, en, curie))
    if cand:
        cand.sort(reverse=True)
        _, s, en, curie = cand[0]
        res = (curie, s == 0 and en >= len(label))   # full-label coverage -> exact
    _ols_cache[key] = res
    response = json.dumps(d, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    _ols_evidence[key] = {
        "query": label,
        "ontology": ont,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "response_sha256": "sha256:"
        + hashlib.sha256(response.encode("utf-8")).hexdigest(),
        "ontology_term": res[0],
        "exact": res[1],
    }
    return res


def ground_label(label: str, ont: str):
    """Ground a label to `ont` (cl|uberon); tries the label and a singular variant, prefers exact."""
    cands, seen = [], set()
    for c in (label, _singularize(label)):
        if c not in seen:
            seen.add(c)
            cands.append(c)
    hits = [h for h in (_ols_once(c, ont) for c in cands) if h[0]]
    if not hits:
        return (None, None)
    hits.sort(key=lambda h: h[1] is True, reverse=True)
    return hits[0]


def ground_term(term: dict, gene_map: dict[str, str | None]) -> dict:
    tt = term.get("term_type")
    lab = term.get("normalized_label") or ""
    if tt == "gene":
        normalized_id = resolve_gene_id(lab, gene_map)
        source_id = resolve_gene_id(term.get("sub_span") or "", gene_map)
        if normalized_id and source_id and normalized_id != source_id:
            raise GeneGroundingConflict(
                f"gene label {lab!r} resolves to {normalized_id}, but source span "
                f"{term.get('sub_span')!r} resolves to {source_id}"
            )
        term["ontology_term"] = normalized_id
        term["exact"] = True if term["ontology_term"] is not None else None
    elif tt in ("celltype", "comparison"):
        term["ontology_term"], term["exact"] = ground_label(lab, "cl")
    elif tt == "tissue":
        term["ontology_term"], term["exact"] = ground_label(lab, "uberon")
    # Preserve the exact source offset within span_literal.
    ss, span = term.get("sub_span"), term.get("_span_literal")
    if ss and span is not None:
        i = span.find(ss)
        term["sub_offset"] = [i, i + len(ss)] if i >= 0 else None
    return term


def ground_claims(claims: list, gene_map_file: Optional[str] = None) -> list:
    """Ground every term in every claim. Loads the offline gene map once."""
    gene_map = load_gene_map(
        gene_map_file=Path(gene_map_file) if gene_map_file else None
    )
    for c in claims:
        span = c.get("span_literal", "")
        for t in c.get("terms", []):
            t["_span_literal"] = span
            ground_term(t, gene_map)
            t.pop("_span_literal", None)
    return claims


def _gene_map_metadata(gene_map_file: Optional[str], organism: str) -> dict:
    path = Path(gene_map_file) if gene_map_file else Path(__file__).parent / "data" / "gmap.txt"
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "provider": "offline-gene-map",
        "organism": organism,
        "sha256": f"sha256:{digest}",
    }


def _document_ols_evidence(document: dict) -> list[dict]:
    evidence: list[dict] = []
    seen: set[tuple[str, str]] = set()
    ontology_by_type = {"celltype": "cl", "comparison": "cl", "tissue": "uberon"}
    for claim in document["claims"]:
        for term in claim["terms"]:
            ontology = ontology_by_type.get(term.get("term_type"))
            if ontology is None:
                continue
            label = term.get("normalized_label") or ""
            for query in dict.fromkeys((label, _singularize(label))):
                key = (query.lower(), ontology)
                if key in seen or key not in _ols_evidence:
                    continue
                seen.add(key)
                evidence.append(copy.deepcopy(_ols_evidence[key]))
    return evidence


def ground_document(
    document: dict,
    *,
    organism: str,
    manuscript_text: Optional[str] = None,
    gene_map_file: Optional[str] = None,
) -> dict:
    """Ground a validated claim document and return a validated onto document."""

    if not organism:
        raise ValueError("organism is required for gene grounding")
    if gene_map_file is None and organism != _PACKAGED_GENE_MAP_ORGANISM:
        raise ValueError(
            "the packaged gene map supports only homo_sapiens; "
            "provide --gene-map for another organism"
        )
    assert_valid_document(document, manuscript_text)
    grounded = copy.deepcopy(document)
    ground_claims(grounded["claims"], gene_map_file=gene_map_file)
    grounded["schema_version"] = ONTO_SCHEMA
    grounded["grounding"] = {
        "genes": _gene_map_metadata(gene_map_file, organism),
        "ontology_service": {
            "provider": "OLS4",
            "endpoint": _OLS_URL,
            "queries": _document_ols_evidence(grounded),
        },
    }
    assert_valid_document(grounded, manuscript_text)
    return grounded
