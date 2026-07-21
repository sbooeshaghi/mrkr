"""Offline gene symbol to Ensembl ID resolution."""

import unicodedata
from pathlib import Path

_GREEK_CHAR_MAP = str.maketrans(
    {
        "Α": "A",
        "α": "A",
        "Β": "B",
        "β": "B",
        "Γ": "G",
        "γ": "G",
    }
)

_DASH_CHAR_MAP = str.maketrans(
    {
        "‐": "-",
        "‑": "-",
        "‒": "-",
        "–": "-",
        "—": "-",
        "−": "-",
    }
)

# Narrow alias table for common human marker synonyms and OCR-style variants.
# Ambiguous protein complexes such as CD3 or HLA-DR are intentionally left
# unmapped.
_GENE_ALIASES = {
    "ADAR-P150": "ADAR",
    "BDCA-2": "CLEC4C",
    "DESMIN": "DES",
    "DNASE13": "DNASE1L3",
    "ECAD": "CDH1",
    "KI67": "MKI67",
    "MIK67": "MKI67",
    "NEPHRIN": "NPHS1",
    "PECAM": "PECAM1",
    "PDGRB": "PDGFRB",
    "RSG10": "RGS10",
    "SCL17A7": "SLC17A7",
    "VISG4": "VSIG4",
}


def _normalize_gene_key(gene_name: str) -> str:
    """Normalize a gene label for lookup."""
    key = unicodedata.normalize("NFKC", gene_name or "")
    key = key.translate(_DASH_CHAR_MAP)
    key = key.translate(_GREEK_CHAR_MAP)
    key = key.upper().strip()
    key = " ".join(key.split())
    return key


def _candidate_gene_keys(gene_name: str) -> list[str]:
    """Generate lookup keys for a reported gene label."""
    base = _normalize_gene_key(gene_name)
    if not base:
        return []

    candidates = []

    def add(key: str) -> None:
        if key and key not in candidates:
            candidates.append(key)

    add(base)

    compact = base.replace(" ", "")
    add(compact)

    if "-" in compact:
        add(compact.replace("-", ""))

    for key in list(candidates):
        alias = _GENE_ALIASES.get(key)
        if alias:
            add(_normalize_gene_key(alias))

    return candidates


def resolve_gene_id(gene_name: str, gene_map: dict[str, str | None]) -> str | None:
    """Resolve a reported gene label to an Ensembl ID using mrkr lookup rules."""
    for key in _candidate_gene_keys(gene_name):
        ensembl_id = gene_map.get(key)
        if ensembl_id:
            return ensembl_id
    return None


def load_gene_map(gene_map_file: Path | None = None) -> dict[str, str | None]:
    """Load a symbol map and leave conflicting symbols unresolved."""
    if gene_map_file is None:
        # Use packaged gmap.txt
        package_dir = Path(__file__).parent
        gene_map_file = package_dir / "data" / "gmap.txt"

    if not gene_map_file.exists():
        raise FileNotFoundError(f"Gene mapping file not found: {gene_map_file}")

    gene_map: dict[str, str | None] = {}
    with open(gene_map_file, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) != 2:
                continue

            gene_name, ensembl_id = parts
            key = _normalize_gene_key(gene_name)
            previous = gene_map.get(key)
            if key in gene_map and previous != ensembl_id:
                gene_map[key] = None
            else:
                gene_map[key] = ensembl_id

    return gene_map
