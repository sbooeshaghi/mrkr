"""Offline tests for mrkr.ground (no network: gene map + validation + singularization)."""

from mrkr.ground import _singularize, ground_term, validate
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
    assert t["direction"] == "positive"             # default filled
    assert t["sub_offset"] == [14, 18]              # index of "CD14" in the span


def test_validate():
    claims = [{
        "span_literal": "Macrophages express CD14.",
        "summary": "The macrophage expresses CD14.",
        "terms": [
            {"sub_span": "Macrophages", "normalized_label": "macrophage", "term_type": "celltype"},
            {"sub_span": "CD14", "normalized_label": "CD14", "term_type": "gene"},
        ],
    }]
    r = validate(claims, paper_text="... Macrophages express CD14. ...")
    assert r["span_ok"] == r["span_total"] == 1
    assert r["sub_ok"] == r["sub_total"] == 2
    assert r["label_ok"] == r["label_total"] == 2

    bad = [{"span_literal": "Macrophages express CD14.", "summary": "no gene here",
            "terms": [{"sub_span": "NOTPRESENT", "normalized_label": "macrophage", "term_type": "celltype"}]}]
    rb = validate(bad, paper_text="Macrophages express CD14.")
    assert rb["sub_ok"] == 0 and rb["label_ok"] == 0


if __name__ == "__main__":
    test_singularize()
    test_gene_grounding_offline()
    test_validate()
    print("ok: all ground tests passed")
