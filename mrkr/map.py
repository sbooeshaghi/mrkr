"""Gene name to Ensembl ID mapping functionality."""

import json
import unicodedata
from pathlib import Path
from typing import Optional

import click


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


def resolve_gene_id(gene_name: str, gene_map: dict[str, str]) -> str | None:
    """Resolve a reported gene label to an Ensembl ID using mrkr lookup rules."""
    for key in _candidate_gene_keys(gene_name):
        ensembl_id = gene_map.get(key)
        if ensembl_id:
            return ensembl_id
    return None


def load_gene_map(gene_map_file: Optional[Path] = None, verbose: bool = False) -> dict[str, str]:
    """
    Load gene name to Ensembl ID mapping.

    Args:
        gene_map_file: Path to custom gene mapping file (tab-separated: gene_name\tensembl_id)
        verbose: Print loading information

    Returns:
        Dictionary mapping uppercase gene names to Ensembl IDs
    """
    if gene_map_file is None:
        # Use packaged gmap.txt
        package_dir = Path(__file__).parent
        gene_map_file = package_dir / "data" / "gmap.txt"

    if not gene_map_file.exists():
        raise FileNotFoundError(f"Gene mapping file not found: {gene_map_file}")

    if verbose:
        click.echo(f"📖 Loading gene mapping from: {gene_map_file}")

    gene_map = {}
    with open(gene_map_file, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) != 2:
                if verbose:
                    click.echo(f"⚠️  Skipping malformed line {line_num}: {line}")
                continue

            gene_name, ensembl_id = parts
            gene_map[_normalize_gene_key(gene_name)] = ensembl_id

    if verbose:
        click.echo(f"✅ Loaded {len(gene_map):,} gene mappings")

    return gene_map


def map_gene_ids(
    input_file: Path,
    gene_map_file: Optional[Path] = None,
    verbose: bool = False,
) -> list[dict]:
    """
    Map gene names to Ensembl IDs in a markers JSON file.

    Args:
        input_file: Path to markers JSON file (output from mrkr extract)
        gene_map_file: Optional custom gene mapping file
        verbose: Print progress information

    Returns:
        List of marker records with feature_id populated where possible
    """
    # Load input markers
    if verbose:
        click.echo(f"📖 Loading markers from: {input_file}")

    with open(input_file, "r", encoding="utf-8") as f:
        markers = json.load(f)

    if not isinstance(markers, list):
        raise ValueError(f"Expected JSON array, got {type(markers).__name__}")

    if verbose:
        click.echo(f"✅ Loaded {len(markers)} marker records")

    # Load gene mapping
    gene_map = load_gene_map(gene_map_file=gene_map_file, verbose=verbose)

    # Map feature_name to feature_id
    if verbose:
        click.echo("\n🔄 Mapping gene names to Ensembl IDs...")

    eligible_count = 0
    mapped_count = 0
    unmapped_genes = set()

    for record in markers:
        feature_name = record.get("feature_name", "")
        if not feature_name:
            continue

        # Only map human genes (gmap.txt is human-only)
        organism = (record.get("organism") or "").strip().lower()
        if organism and organism != "homo_sapiens":
            continue

        eligible_count += 1

        ensembl_id = resolve_gene_id(feature_name, gene_map)

        if ensembl_id:
            record["feature_id"] = ensembl_id
            mapped_count += 1
        else:
            unmapped_genes.add(feature_name)

    if verbose:
        if eligible_count == 0:
            click.echo("✅ No eligible human marker genes found; nothing to map")
        else:
            pct = mapped_count / eligible_count * 100
            click.echo(f"✅ Mapped {mapped_count}/{eligible_count} genes ({pct:.1f}%)")
        if unmapped_genes:
            click.echo(f"⚠️  {len(unmapped_genes)} unique genes could not be mapped:")
            for gene in sorted(unmapped_genes)[:10]:
                click.echo(f"   - {gene}")
            if len(unmapped_genes) > 10:
                click.echo(f"   ... and {len(unmapped_genes) - 10} more")

    return markers
