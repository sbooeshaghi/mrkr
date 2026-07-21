"""Offline tests for mrkr.ground."""

import urllib.error

import pytest

from mrkr.ground import (
    GeneGroundingConflict,
    GroundingServiceError,
    _ols_cache,
    _ols_once,
    _singularize,
    ground_term,
)
from mrkr.map import load_gene_map


def test_singularize():
    assert _singularize("oocytes") == "oocyte"
    assert _singularize("granulosa cells") == "granulosa cell"
    assert _singularize("bodies") == "body"
    assert _singularize("mass") == "mass"          # -ss not stripped
    assert _singularize("T-cells") == "T-cell"


def test_gene_grounding_offline():
    gene_map = load_gene_map()
    t = {"term_type": "gene", "normalized_label": "CD14", "sub_span": "CD14",
         "_span_literal": "cells express CD14 strongly"}
    ground_term(t, gene_map)
    assert t["ontology_term"] and t["ontology_term"].startswith("ENSG")
    assert t["exact"] is True
    assert "direction" not in t
    assert t["sub_offset"] == [14, 18]              # index of "CD14" in the span


def test_unresolved_gene_has_unknown_exactness():
    term = {"term_type": "gene", "normalized_label": "NOT_A_REAL_GENE"}

    ground_term(term, {})

    assert term["ontology_term"] is None
    assert term["exact"] is None


def test_gene_grounding_rejects_source_normalization_conflict():
    term = {
        "term_type": "gene",
        "normalized_label": "LYZ",
        "sub_span": "CD14",
    }
    gene_map = {"LYZ": "ENSG_LYZ", "CD14": "ENSG_CD14"}

    with pytest.raises(GeneGroundingConflict):
        ground_term(term, gene_map)


def test_ols_failure_is_not_reported_as_unresolved(monkeypatch):
    _ols_cache.clear()

    def fail(*args, **kwargs):
        raise urllib.error.URLError("offline")

    monkeypatch.setattr("urllib.request.urlopen", fail)
    with pytest.raises(GroundingServiceError, match="OLS request failed"):
        _ols_once("macrophage", "cl")


if __name__ == "__main__":
    test_singularize()
    test_gene_grounding_offline()
    print("ok: all ground tests passed")
