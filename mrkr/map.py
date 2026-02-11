"""Gene name to Ensembl ID mapping functionality."""

import json
from pathlib import Path
from typing import Optional

import click


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
            # Store in uppercase for case-insensitive matching
            gene_map[gene_name.upper()] = ensembl_id

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

    mapped_count = 0
    unmapped_genes = set()

    for record in markers:
        feature_name = record.get("feature_name", "")
        if not feature_name:
            continue

        # Look up in gene map (already uppercase)
        ensembl_id = gene_map.get(feature_name.upper())

        if ensembl_id:
            record["feature_id"] = ensembl_id
            mapped_count += 1
        else:
            unmapped_genes.add(feature_name)

    if verbose:
        click.echo(f"✅ Mapped {mapped_count}/{len(markers)} genes ({mapped_count/len(markers)*100:.1f}%)")
        if unmapped_genes:
            click.echo(f"⚠️  {len(unmapped_genes)} unique genes could not be mapped:")
            for gene in sorted(unmapped_genes)[:10]:
                click.echo(f"   - {gene}")
            if len(unmapped_genes) > 10:
                click.echo(f"   ... and {len(unmapped_genes) - 10} more")

    return markers
