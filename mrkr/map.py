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
    "C-KIT": "KIT",
    "CKIT": "KIT",
    "DESMIN": "DES",
    "DNASE13": "DNASE1L3",
    "ECAD": "CDH1",
    "INTERFERON-GAMMA": "IFNG",
    "INTERLEUKIN-1 RECEPTOR TYPE 1": "IL1R1",
    "INTERLEUKIN-13": "IL13",
    "INTERLEUKIN-17A": "IL17A",
    "INTERLEUKIN-2": "IL2",
    "IGFB3": "IGFBP3",
    "KI67": "MKI67",
    "LIPOPROTEIN LIPASE": "LPL",
    "MIK67": "MKI67",
    "NG2": "CSPG4",
    "NEPHRIN": "NPHS1",
    "NKP44": "NCR2",
    "NOV": "CCN3",
    "PECAM": "PECAM1",
    "PDGRB": "PDGFRB",
    "PERFORIN": "PRF1",
    "RSG10": "RGS10",
    "SCL17A7": "SLC17A7",
    "T-BET": "TBX21",
    "TBET": "TBX21",
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


def _load_rows(gene_map_file: Path) -> dict[str, str | None]:
    """Load one two-column map and leave conflicting labels unresolved."""

    gene_map: dict[str, str | None] = {}
    with gene_map_file.open("r", encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            if not line or line.startswith("#"):
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


def load_gene_map(
    gene_map_file: Path | None = None,
    canonical_gene_map_file: Path | None = None,
) -> dict[str, str | None]:
    """Load aliases from gmap.txt, preferring authoritative gene symbols."""

    if gene_map_file is None:
        package_dir = Path(__file__).parent
        gene_map_file = package_dir / "data" / "gmap.txt"
        canonical_gene_map_file = package_dir / "data" / "gmap_canonical.txt"

    if not gene_map_file.exists():
        raise FileNotFoundError(f"Gene mapping file not found: {gene_map_file}")

    gene_map = _load_rows(gene_map_file)
    if canonical_gene_map_file is not None:
        if not canonical_gene_map_file.exists():
            raise FileNotFoundError(
                f"Canonical gene mapping file not found: {canonical_gene_map_file}"
            )
        canonical = _load_rows(canonical_gene_map_file)
        for key, ensembl_id in canonical.items():
            if ensembl_id is not None:
                gene_map[key] = ensembl_id

    return gene_map
