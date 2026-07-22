"""Tests for the canonical mrkr claim and onto document contract."""

import copy

import pytest

from mrkr.claims import (
    CLAIMS_SCHEMA,
    ONTO_SCHEMA,
    ClaimValidationError,
    assert_valid_document,
    make_claim_document,
    prepare_raw_claims,
    validate_document,
)
from mrkr.ground import ground_document

MANUSCRIPT = "Macrophages express CD14 in human blood."
RAW_CLAIMS = [
    {
        "span_literal": MANUSCRIPT,
        "summary": "In Homo sapiens, macrophage expresses CD14 in blood.",
        "terms": [
            {
                "sub_span": "Macrophages",
                "normalized_label": "macrophage",
                "term_type": "celltype",
            },
            {
                "sub_span": "CD14",
                "normalized_label": "CD14",
                "term_type": "gene",
                "direction": "positive",
            },
            {
                "sub_span": "blood",
                "normalized_label": "blood",
                "term_type": "tissue",
            },
            {
                "sub_span": "human",
                "normalized_label": "Homo sapiens",
                "term_type": "organism",
            },
        ],
    }
]


def test_make_and_validate_claim_document():
    document = make_claim_document(
        source_id="paper.md", manuscript_text=MANUSCRIPT, raw_claims=RAW_CLAIMS
    )
    assert document["schema_version"] == CLAIMS_SCHEMA
    assert document["claims"][0]["span_offset"] == [0, len(MANUSCRIPT)]
    assert document["claims"][0]["terms"][1]["sub_offset"] == [20, 24]
    assert_valid_document(document, MANUSCRIPT)


def test_validation_reports_cardinality_and_alignment_errors():
    document = make_claim_document(
        source_id="paper.md", manuscript_text=MANUSCRIPT, raw_claims=RAW_CLAIMS
    )
    broken = copy.deepcopy(document)
    broken["claims"][0]["terms"] = broken["claims"][0]["terms"][1:]
    broken["claims"][0]["terms"][0]["sub_offset"] = [0, 4]
    report = validate_document(broken, MANUSCRIPT)
    codes = {error["code"] for error in report["errors"]}
    assert report["valid"] is False
    assert "claim.celltype_count" in codes
    assert "term.sub_offset" in codes
    with pytest.raises(ClaimValidationError):
        assert_valid_document(broken, MANUSCRIPT)


def test_validation_rejects_wrong_source_hash():
    document = make_claim_document(
        source_id="paper.md", manuscript_text=MANUSCRIPT, raw_claims=RAW_CLAIMS
    )
    report = validate_document(document, MANUSCRIPT + " Changed.")
    assert any(error["code"] == "source.hash_mismatch" for error in report["errors"])


def test_repeated_explicit_evidence_uses_a_valid_occurrence():
    repeated = f"{MANUSCRIPT}\n{MANUSCRIPT}"
    document = make_claim_document(
        source_id="paper.md", manuscript_text=repeated, raw_claims=RAW_CLAIMS
    )
    assert_valid_document(document, repeated)


def test_repeated_evidence_with_implicit_target_is_ambiguous():
    repeated = f"{MANUSCRIPT}\n{MANUSCRIPT}"
    implicit = copy.deepcopy(RAW_CLAIMS)
    implicit[0]["terms"][0]["sub_span"] = None
    document = make_claim_document(
        source_id="paper.md", manuscript_text=repeated, raw_claims=implicit
    )
    report = validate_document(document, repeated)
    assert any(error["code"] == "claim.span_offset" for error in report["errors"])


def test_same_span_and_target_become_one_marker_panel():
    manuscript = "Macrophages express CD14 and LYZ in human blood."
    first = copy.deepcopy(RAW_CLAIMS[0])
    first["span_literal"] = manuscript
    second = copy.deepcopy(first)
    second["summary"] = "macrophage expresses LYZ in blood."
    second["terms"][1] = {
        "sub_span": "LYZ",
        "normalized_label": "LYZ",
        "term_type": "gene",
        "direction": "positive",
    }
    document = make_claim_document(
        source_id="paper.md",
        manuscript_text=manuscript,
        raw_claims=[first, second],
    )
    assert len(document["claims"]) == 1
    genes = [
        term["normalized_label"]
        for term in document["claims"][0]["terms"]
        if term["term_type"] == "gene"
    ]
    assert genes == ["CD14", "LYZ"]
    assert_valid_document(document, manuscript)


def test_validation_rejects_implicit_marker_gene():
    raw = copy.deepcopy(RAW_CLAIMS)
    raw[0]["terms"][1]["sub_span"] = None
    document = make_claim_document(
        source_id="paper.md", manuscript_text=MANUSCRIPT, raw_claims=raw
    )

    report = validate_document(document, MANUSCRIPT)

    assert any(error["code"] == "term.gene_explicit" for error in report["errors"])


def test_prepare_raw_claims_reanchors_omitted_parenthetical():
    manuscript = (
        "Macrophages (Figure 2) express CD14 and LYZ in human blood."
    )
    raw = copy.deepcopy(RAW_CLAIMS)
    raw[0]["span_literal"] = "Macrophages express CD14 and LYZ in human blood."
    raw[0]["terms"].insert(
        2,
        {
            "sub_span": "LYZ",
            "normalized_label": "LYZ",
            "term_type": "gene",
            "direction": "positive",
        },
    )

    prepared, report = prepare_raw_claims(manuscript, raw)

    assert prepared[0]["span_literal"] == manuscript
    assert report["reanchored_spans"] == 1


def test_prepare_raw_claims_excludes_records_without_a_gene():
    raw = copy.deepcopy(RAW_CLAIMS)
    raw[0]["terms"] = [raw[0]["terms"][0]]

    prepared, report = prepare_raw_claims(MANUSCRIPT, raw)

    assert prepared == []
    assert report["excluded_claims"] == [
        {"raw_index": 0, "reason": "no_marker_gene"}
    ]


def test_prepare_raw_claims_excludes_claims_without_an_explicit_gene():
    raw = copy.deepcopy(RAW_CLAIMS)
    raw[0]["terms"][1]["sub_span"] = None

    prepared, report = prepare_raw_claims(MANUSCRIPT, raw)

    assert prepared == []
    assert report["excluded_claims"] == [
        {"raw_index": 0, "reason": "no_explicit_marker_gene"}
    ]
    assert report["excluded_terms"] == [
        {
            "raw_index": 0,
            "normalized_label": "CD14",
            "reason": "implicit_marker_gene",
        }
    ]


def test_prepare_raw_claims_expands_gene_family_shorthand():
    manuscript = "Cells are enriched for FOXP1/2 regulons."
    raw = copy.deepcopy(RAW_CLAIMS)
    raw[0]["span_literal"] = manuscript
    raw[0]["terms"][0]["sub_span"] = "Cells"
    raw[0]["terms"][1] = {
        "sub_span": "FOXP2",
        "normalized_label": "FOXP2",
        "term_type": "gene",
        "direction": "positive",
    }

    prepared, report = prepare_raw_claims(manuscript, raw)

    assert prepared[0]["terms"][1]["sub_span"] == "FOXP1/2"
    assert report["expanded_gene_shorthand"] == 1


def test_prepare_raw_claims_reconstructs_formatted_marker_clause():
    manuscript = (
        "Flow cytometry identified TBET^(hi)EOMES^(int)PERFORIN^(hi)CD16+ "
        "mNK cells in tissue."
    )
    raw = [
        {
            "span_literal": "TBET+CD16+ mNK cells",
            "summary": (
                "In Homo sapiens, mature natural killer cells express TBET and CD16."
            ),
            "terms": [
                {
                    "sub_span": "mNK cells",
                    "normalized_label": "mature natural killer cell",
                    "term_type": "celltype",
                },
                {
                    "sub_span": "TBET",
                    "normalized_label": "TBET",
                    "term_type": "gene",
                    "direction": "positive",
                },
                {
                    "sub_span": "CD16",
                    "normalized_label": "CD16",
                    "term_type": "gene",
                    "direction": "positive",
                },
                {
                    "sub_span": None,
                    "normalized_label": "Homo sapiens",
                    "term_type": "organism",
                },
            ],
        }
    ]

    prepared, report = prepare_raw_claims(manuscript, raw)
    document = make_claim_document(
        source_id="paper.md", manuscript_text=manuscript, raw_claims=prepared
    )

    assert prepared[0]["span_literal"] == (
        "TBET^(hi)EOMES^(int)PERFORIN^(hi)CD16+ mNK cells"
    )
    assert report["reconstructed_spans"] == 1
    assert_valid_document(document, manuscript)


def test_prepare_raw_claims_reconstruction_keeps_context_terms():
    manuscript = "WAT cDC1s express DPP4 and are CD1C^(-)."
    raw = [
        {
            "span_literal": "WAT cDC1s express DPP4 and are CD1C(-).",
            "summary": (
                "In Homo sapiens, white adipose tissue conventional type 1 dendritic cell "
                "expresses DPP4, not CD1C."
            ),
            "terms": [
                {
                    "sub_span": "cDC1s",
                    "normalized_label": "conventional type 1 dendritic cell",
                    "term_type": "celltype",
                },
                {
                    "sub_span": "DPP4",
                    "normalized_label": "DPP4",
                    "term_type": "gene",
                    "direction": "positive",
                },
                {
                    "sub_span": "CD1C",
                    "normalized_label": "CD1C",
                    "term_type": "gene",
                    "direction": "negative",
                },
                {
                    "sub_span": "WAT",
                    "normalized_label": "white adipose tissue",
                    "term_type": "tissue",
                },
                {
                    "sub_span": None,
                    "normalized_label": "Homo sapiens",
                    "term_type": "organism",
                },
            ],
        }
    ]

    prepared, report = prepare_raw_claims(manuscript, raw)
    document = make_claim_document(
        source_id="paper.md", manuscript_text=manuscript, raw_claims=prepared
    )

    assert prepared[0]["span_literal"] == "WAT cDC1s express DPP4 and are CD1C"
    assert report["reconstructed_spans"] == 1
    assert_valid_document(document, manuscript)


def test_ground_document_keeps_unresolved_terms_explicit(monkeypatch):
    document = make_claim_document(
        source_id="paper.md", manuscript_text=MANUSCRIPT, raw_claims=RAW_CLAIMS
    )

    def fake_ground_label(label, ontology):
        return {
            ("macrophage", "cl"): ("CL:0000235", True),
            ("blood", "uberon"): (None, None),
        }[(label, ontology)]

    monkeypatch.setattr("mrkr.ground.ground_label", fake_ground_label)
    monkeypatch.setattr(
        "mrkr.ground._document_ols_evidence",
        lambda _document: [
            {
                "query": "macrophage",
                "ontology": "cl",
                "retrieved_at": "2026-01-01T00:00:00+00:00",
                "response_sha256": "sha256:test",
                "ontology_term": "CL:0000235",
                "exact": True,
            },
            {
                "query": "blood",
                "ontology": "uberon",
                "retrieved_at": "2026-01-01T00:00:00+00:00",
                "response_sha256": "sha256:test",
                "ontology_term": None,
                "exact": None,
            },
        ],
    )
    grounded = ground_document(
        document, organism="homo_sapiens", manuscript_text=MANUSCRIPT
    )
    assert grounded["schema_version"] == ONTO_SCHEMA
    terms = grounded["claims"][0]["terms"]
    assert terms[0]["ontology_term"] == "CL:0000235"
    assert terms[1]["ontology_term"].startswith("ENSG")
    assert terms[2]["ontology_term"] is None
    assert terms[2]["exact"] is None
    assert_valid_document(grounded, MANUSCRIPT)


def test_onto_document_requires_grounding_fields_and_metadata():
    document = make_claim_document(
        source_id="paper.md", manuscript_text=MANUSCRIPT, raw_claims=RAW_CLAIMS
    )
    document["schema_version"] = ONTO_SCHEMA

    report = validate_document(document, MANUSCRIPT)

    codes = {error["code"] for error in report["errors"]}
    assert "grounding.missing" in codes
    assert "term.grounding_state" in codes


def test_onto_document_rejects_empty_grounding_metadata():
    document = make_claim_document(
        source_id="paper.md", manuscript_text=MANUSCRIPT, raw_claims=RAW_CLAIMS
    )
    for term in document["claims"][0]["terms"]:
        term["ontology_term"] = None
        term["exact"] = None
    document["schema_version"] = ONTO_SCHEMA
    document["grounding"] = {}

    report = validate_document(document, MANUSCRIPT)
    codes = {error["code"] for error in report["errors"]}

    assert "grounding.genes" in codes
    assert "grounding.service" in codes


def test_validation_reports_errors_and_grounding_warnings(monkeypatch):
    document = make_claim_document(
        source_id="paper.md", manuscript_text=MANUSCRIPT, raw_claims=RAW_CLAIMS
    )

    def fake_ground_label(label, ontology):
        return {
            ("macrophage", "cl"): ("CL:0000235", False),
            ("blood", "uberon"): (None, None),
        }[(label, ontology)]

    monkeypatch.setattr("mrkr.ground.ground_label", fake_ground_label)
    monkeypatch.setattr(
        "mrkr.ground._document_ols_evidence",
        lambda _document: [
            {
                "query": "macrophage",
                "ontology": "cl",
                "retrieved_at": "2026-01-01T00:00:00+00:00",
                "response_sha256": "sha256:test",
                "ontology_term": "CL:0000235",
                "exact": False,
            },
            {
                "query": "blood",
                "ontology": "uberon",
                "retrieved_at": "2026-01-01T00:00:00+00:00",
                "response_sha256": "sha256:test",
                "ontology_term": None,
                "exact": None,
            },
        ],
    )
    grounded = ground_document(
        document, organism="homo_sapiens", manuscript_text=MANUSCRIPT
    )

    report = validate_document(grounded, MANUSCRIPT)

    assert report["valid"] is True
    assert report["n_errors"] == 0
    assert {warning["code"] for warning in report["warnings"]} == {
        "term.grounding_coarse",
        "term.grounding_unresolved",
    }


def test_validation_warns_when_nonorganism_label_is_missing_from_summary():
    document = make_claim_document(
        source_id="paper.md", manuscript_text=MANUSCRIPT, raw_claims=RAW_CLAIMS
    )
    document["claims"][0]["summary"] = "In Homo sapiens, a marker statement."

    report = validate_document(document, MANUSCRIPT)

    assert report["valid"] is True
    assert any(
        warning["code"] == "term.label_missing_summary"
        for warning in report["warnings"]
    )


def test_validation_requires_organism_label_in_summary():
    document = make_claim_document(
        source_id="paper.md", manuscript_text=MANUSCRIPT, raw_claims=RAW_CLAIMS
    )
    document["claims"][0]["summary"] = "Macrophage expresses CD14 in blood."

    report = validate_document(document, MANUSCRIPT)

    assert report["valid"] is False
    assert any(
        error["code"] == "term.label_missing_summary"
        for error in report["errors"]
    )
