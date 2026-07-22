"""Canonical organism names used by extraction and grounding."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Organism:
    """One supported organism and its NCBI Taxonomy identifier."""

    key: str
    label: str
    ontology_term: str


_ORGANISMS = (
    Organism("homo_sapiens", "Homo sapiens", "NCBITaxon:9606"),
    Organism("mus_musculus", "Mus musculus", "NCBITaxon:10090"),
    Organism("macaca_mulatta", "Macaca mulatta", "NCBITaxon:9544"),
    Organism("danio_rerio", "Danio rerio", "NCBITaxon:7955"),
    Organism("caenorhabditis_elegans", "Caenorhabditis elegans", "NCBITaxon:6239"),
    Organism("drosophila_melanogaster", "Drosophila melanogaster", "NCBITaxon:7227"),
    Organism("rattus_norvegicus", "Rattus norvegicus", "NCBITaxon:10116"),
    Organism("saccharomyces_cerevisiae", "Saccharomyces cerevisiae", "NCBITaxon:4932"),
    Organism("arabidopsis_thaliana", "Arabidopsis thaliana", "NCBITaxon:3702"),
)


def _key(value: str) -> str:
    return "_".join(value.strip().lower().replace("-", " ").split())


_BY_KEY = {
    alias: organism
    for organism in _ORGANISMS
    for alias in (organism.key, _key(organism.label))
}
_BY_KEY.update(
    {
        "human": _BY_KEY["homo_sapiens"],
        "mouse": _BY_KEY["mus_musculus"],
        "macaque": _BY_KEY["macaca_mulatta"],
        "zebrafish": _BY_KEY["danio_rerio"],
        "worm": _BY_KEY["caenorhabditis_elegans"],
        "fruit_fly": _BY_KEY["drosophila_melanogaster"],
        "rat": _BY_KEY["rattus_norvegicus"],
        "yeast": _BY_KEY["saccharomyces_cerevisiae"],
    }
)


def get_organism(value: str) -> Organism:
    """Return the canonical record for a supported organism name or alias."""

    organism = _BY_KEY.get(_key(value))
    if organism is None:
        supported = ", ".join(item.key for item in _ORGANISMS)
        raise ValueError(f"unsupported organism {value!r}; choose one of: {supported}")
    return organism
