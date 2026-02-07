"""CLI interface for mrkr."""

import json
import sys
from pathlib import Path
from typing import Optional, Tuple

import click

from .config import config
from .extract import extract_markers
from .map import map_gene_ids


def parse_file_list(
    ctx: click.Context, param: click.Parameter, value: Optional[Tuple[str, ...]]
) -> Optional[list[Path]]:
    r"""Parse file arguments with support for glob patterns.

    Supports:
    - Multiple -f flags: -f file1.tif -f file2.tif -f file3.tif
    - Glob patterns: -f "*.tif" or -f "fig*.tif"
    - Files with spaces: -f "My File.tif"
    """
    if not value:
        return None

    import glob as glob_module

    paths = []
    for path_str in value:
        # Check if this is a glob pattern
        if any(char in path_str for char in ['*', '?', '[', ']']):
            # Expand glob pattern
            expanded = glob_module.glob(path_str)
            if not expanded:
                raise click.BadParameter(f"No files match pattern: {path_str}")
            for expanded_path in sorted(expanded):
                path = Path(expanded_path)
                if path.exists():
                    paths.append(path)
        else:
            # Regular path
            path = Path(path_str)
            if not path.exists():
                raise click.BadParameter(f"File not found: {path}")
            paths.append(path)

    return paths if paths else None


@click.group()
@click.version_option(version="0.2.0", prog_name="mrkr")
def cli():
    """
    mrkr - Cell type marker gene extraction and mapping tool.

    Commands:
        extract - Extract marker genes from manuscripts, figures, and DEG tables
        map     - Map gene names to Ensembl IDs
    """
    pass


@cli.command()
@click.option(
    "--manuscript",
    "-m",
    type=click.Path(exists=True, path_type=Path),
    help="Manuscript markdown file",
)
@click.option(
    "--figures",
    "-f",
    "figures",
    type=str,
    multiple=True,
    callback=parse_file_list,
    help="Figure files. Use: -f file1.png -f file2.png OR -f '*.png' OR -f 'My File.png'",
)
@click.option(
    "--deg",
    "-d",
    "deg",
    type=str,
    multiple=True,
    callback=parse_file_list,
    help="DEG files. Use: -d file1.xlsx -d file2.xlsx OR -d '*.xlsx' OR -d 'My File.xlsx'",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(path_type=Path),
    help="Output JSON file",
)
@click.option(
    "--species",
    "-s",
    help="Species Latin name (e.g., homo_sapiens). If not provided, inferred from manuscript",
)
@click.option(
    "--metrics",
    type=click.Path(path_type=Path),
    help="Output file for LLM metrics (token usage, costs, timing)",
)
@click.option(
    "--verbose", "-v", is_flag=True, help="Verbose output showing progress and token usage"
)
def extract(manuscript, figures, deg, output, species, metrics, verbose):
    """
    Extract cell type marker genes from manuscripts, figures, and/or DEG tables.

    \b
    Examples:
        # Extract from manuscript only
        mrkr extract -m manuscript.md -o markers.json

        # Extract from manuscript + figures (multiple ways to specify files)
        mrkr extract -m manuscript.md -f fig1.png fig2.png -o markers.json
        mrkr extract -m manuscript.md -f "*.png" -o markers.json
        mrkr extract -m manuscript.md -f "*_fig*.tif" -o markers.json

        # Extract from manuscript + DEG tables (uses DEG cell type names)
        mrkr extract -m manuscript.md -d deg_table.xlsx -o markers.json

        # Extract from DEG table only
        mrkr extract -d deg_table.xlsx -s homo_sapiens -o markers.json

        # Full extraction with all sources
        mrkr extract -m manuscript.md -f fig1.png -f fig2.png -d deg.xlsx -o markers.json

        # Multiple DEG files
        mrkr extract -d deg1.xlsx deg2.csv -s homo_sapiens -o markers.json
        mrkr extract -d "Supplementary Table S*.xlsx" -s homo_sapiens -o markers.json

    \b
    Supported formats:
        - Manuscripts: .md, .txt (markdown or plain text)
        - Figures: .png, .jpg, .jpeg
        - DEG tables: .xlsx, .csv, .tsv
    """
    # Validate inputs
    if not manuscript and not figures and not deg:
        click.echo(
            "❌ Error: At least one of --manuscript, --figures, or --deg must be provided",
            err=True,
        )
        raise click.Abort()

    # Convert tuples to lists
    figure_paths = list(figures) if figures else None
    deg_paths = list(deg) if deg else None

    # Validate configuration
    try:
        config.validate()
    except ValueError as e:
        click.echo(f"❌ Configuration error: {e}", err=True)
        raise click.Abort()

    # Print header
    if not verbose:
        click.echo("\n🧬 mrkr extract - Marker gene extraction")
        click.echo("=" * 50)

    # Run extraction
    try:
        results = extract_markers(
            manuscript_path=manuscript,
            figure_paths=figure_paths,
            deg_paths=deg_paths,
            species=species,
            verbose=verbose,
            metrics_path=metrics,
        )

        # Save results
        output.write_text(json.dumps(results, indent=2), encoding="utf-8")

        if not verbose:
            click.echo(f"\n✅ Extracted {len(results)} marker gene associations")
            click.echo(f"💾 Saved to: {output}")
        else:
            click.echo(f"\n💾 Saved to: {output}")

    except Exception as e:
        click.echo(f"\n❌ Error during extraction: {e}", err=True)
        if verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)


@cli.command()
@click.argument("input", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(path_type=Path),
    help="Output JSON file with mapped gene IDs",
)
@click.option(
    "--gene-map",
    type=click.Path(exists=True, path_type=Path),
    help="Gene mapping file (default: uses packaged gmap.txt)",
)
@click.option(
    "--verbose", "-v", is_flag=True, help="Verbose output"
)
def map(input, output, gene_map, verbose):
    """
    Map gene names to Ensembl IDs.

    Takes a markers JSON file and maps feature_label/feature_name to feature_id
    using the gene mapping file.

    \b
    Examples:
        # Map using default gene mapping
        mrkr map markers.json -o markers_mapped.json

        # Map using custom gene mapping file
        mrkr map markers.json -o markers_mapped.json --gene-map custom_gmap.txt

        # Verbose mode
        mrkr map markers.json -o markers_mapped.json -v

    The gene mapping file should be tab-separated with format:
        gene_name    ensembl_id
        A1BG         ENSG00000121410
        ...
    """
    # Print header
    if not verbose:
        click.echo("\n🧬 mrkr map - Gene ID mapping")
        click.echo("=" * 50)

    # Run mapping
    try:
        results = map_gene_ids(
            input_file=input,
            gene_map_file=gene_map,
            verbose=verbose,
        )

        # Save results
        output.write_text(json.dumps(results, indent=2), encoding="utf-8")

        # Count mapped vs unmapped
        mapped = sum(1 for r in results if r.get("feature_id"))
        total = len(results)

        if not verbose:
            click.echo(f"\n✅ Mapped {mapped}/{total} genes ({mapped/total*100:.1f}%)")
            click.echo(f"💾 Saved to: {output}")
        else:
            click.echo(f"\n💾 Saved to: {output}")

    except Exception as e:
        click.echo(f"\n❌ Error during mapping: {e}", err=True)
        if verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    cli()
