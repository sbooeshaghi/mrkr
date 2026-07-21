"""Tests for offline gene identifier resolution."""

from mrkr.map import load_gene_map, resolve_gene_id


def test_alias_and_unicode_variants_resolve(tmp_path):
    gene_map_path = tmp_path / "gmap.txt"
    gene_map_path.write_text(
        "MKI67\tENSG_MKI67\n"
        "PECAM1\tENSG_PECAM1\n"
        "PDGFRA\tENSG_PDGFRA\n"
        "IFNG\tENSG_IFNG\n"
        "LAG3\tENSG_LAG3\n",
        encoding="utf-8",
    )
    gene_map = load_gene_map(gene_map_path)

    assert resolve_gene_id("KI67", gene_map) == "ENSG_MKI67"
    assert resolve_gene_id("PECAM-1", gene_map) == "ENSG_PECAM1"
    assert resolve_gene_id("PECAM", gene_map) == "ENSG_PECAM1"
    assert resolve_gene_id("PDGFRΑ", gene_map) == "ENSG_PDGFRA"
    assert resolve_gene_id("IFN-Γ", gene_map) == "ENSG_IFNG"
    assert resolve_gene_id("LAG-3", gene_map) == "ENSG_LAG3"
    assert resolve_gene_id("CD3", gene_map) is None


def test_conflicting_gene_symbols_are_unresolved(tmp_path):
    gene_map_path = tmp_path / "gmap.txt"
    gene_map_path.write_text(
        "RBP1\tENSG_ONE\nRBP1\tENSG_TWO\n",
        encoding="utf-8",
    )

    gene_map = load_gene_map(gene_map_path)

    assert resolve_gene_id("RBP1", gene_map) is None
