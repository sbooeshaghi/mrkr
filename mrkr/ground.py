"""Deterministic grounding for mrkr claim objects.

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

import json
import ssl
import urllib.parse
import urllib.request
from typing import Optional

from .map import load_gene_map, resolve_gene_id

_CTX = ssl.create_default_context()
_CTX.check_hostname = False
_CTX.verify_mode = ssl.CERT_NONE
_DELIM = " ,.;:!?\t\n()[]{}\"'/\\-"
_OLS_URL = "https://www.ebi.ac.uk/ols4/api/v2/tag_text"
_ols_cache: dict = {}


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
        with urllib.request.urlopen(req, timeout=30, context=_CTX) as r:
            d = json.load(r)
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
    except Exception:
        pass
    _ols_cache[key] = res
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


def ground_term(term: dict, gene_map: dict) -> dict:
    tt = term.get("term_type")
    lab = term.get("normalized_label") or ""
    if tt == "gene":
        term["ontology_term"] = resolve_gene_id(lab, gene_map)
        term["exact"] = term["ontology_term"] is not None
        if not term.get("direction"):
            term["direction"] = "positive"
    elif tt in ("celltype", "comparison"):
        term["ontology_term"], term["exact"] = ground_label(lab, "cl")
    elif tt == "tissue":
        term["ontology_term"], term["exact"] = ground_label(lab, "uberon")
    # sub_offset within span_literal (deterministic, provenance)
    ss, span = term.get("sub_span"), term.get("_span_literal")
    if ss and span is not None:
        i = span.find(ss)
        term["sub_offset"] = [i, i + len(ss)] if i >= 0 else None
    return term


def ground_claims(claims: list, gene_map_file: Optional[str] = None) -> list:
    """Ground every term in every claim. Loads the offline gene map once."""
    gene_map = load_gene_map(gene_map_file=gene_map_file)
    for c in claims:
        span = c.get("span_literal", "")
        for t in c.get("terms", []):
            t["_span_literal"] = span
            ground_term(t, gene_map)
            t.pop("_span_literal", None)
    return claims


def validate(claims: list, paper_text: Optional[str] = None) -> dict:
    """Provenance checks: span_literal in paper; sub_span in span_literal; normalized_label in summary."""
    r = {"span_ok": 0, "span_total": 0, "sub_ok": 0, "sub_total": 0, "label_ok": 0, "label_total": 0}
    for c in claims:
        span, summ = c.get("span_literal", ""), c.get("summary", "")
        summ_lc = summ.lower()
        if paper_text is not None:
            r["span_total"] += 1
            r["span_ok"] += span in paper_text            # exact: offsets depend on it
        for t in c.get("terms", []):
            ss, nl = t.get("sub_span"), t.get("normalized_label", "")
            if ss is not None:
                r["sub_total"] += 1
                r["sub_ok"] += ss in span                 # exact: sub_offset depends on it
            r["label_total"] += 1
            # case-insensitive: the summary is a canonical rewrite; sentence-initial
            # capitals / plurals ("Macrophages" vs "macrophage") are not real mismatches.
            r["label_ok"] += nl.lower() in summ_lc
    return r
