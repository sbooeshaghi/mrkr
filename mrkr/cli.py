"""CLI interface for mrkr."""

import json
import sys
from pathlib import Path
from typing import Optional, Tuple

import click

from .config import config
from .extract import extract_markers
from .generate import generate_celltypes_to_genes, generate_genes_to_celltypes
from .map import map_gene_ids
from .verify import verify_extractions


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
    "--cell-types",
    "-c",
    type=click.Path(exists=True, path_type=Path),
    help="JSON file with known cell types (e.g., evidence_deg/extracted.json). "
         "Extracts unique group_name values to guide LLM cell type matching.",
)
@click.option(
    "--verify/--no-verify", default=True, help="Verify extractions against manuscript text (default: enabled)"
)
@click.option(
    "--verbose", "-v", is_flag=True, help="Verbose output showing progress and token usage"
)
def extract(manuscript, figures, deg, output, species, metrics, cell_types, verify, verbose):
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

    # Load known cell types from JSON file if provided
    # Builds a dict of {data_id: sorted list of group_names} for per-source context
    known_cell_types_list = None
    if cell_types:
        ct_data = json.loads(cell_types.read_text(encoding="utf-8"))
        if isinstance(ct_data, list):
            from collections import defaultdict
            by_source = defaultdict(set)
            for r in ct_data:
                gn = r.get("group_name")
                did = r.get("data_id", "unknown")
                if gn:
                    by_source[did].add(gn)
            known_cell_types_list = {
                did: sorted(gns) for did, gns in sorted(by_source.items())
            }
        if verbose and known_cell_types_list:
            total = len(set(gn for gns in known_cell_types_list.values() for gn in gns))
            click.echo(f"📋 Loaded {total} known cell types from {len(known_cell_types_list)} sources in {cell_types.name}")

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

    # Build command string for metrics reproducibility
    command_str = " ".join(sys.argv)

    # Run extraction
    try:
        results = extract_markers(
            manuscript_path=manuscript,
            figure_paths=figure_paths,
            deg_paths=deg_paths,
            species=species,
            verbose=verbose,
            metrics_path=metrics,
            verify=verify,
            known_cell_types=known_cell_types_list,
            command=command_str,
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


@cli.command()
@click.argument("input", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--manuscript",
    "-m",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Manuscript text file to verify against",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    help="Output JSON file with verification results (default: overwrites input)",
)
@click.option(
    "--verbose", "-v", is_flag=True, help="Verbose output"
)
def verify(input, manuscript, output, verbose):
    """
    Verify extractions against manuscript text.

    Checks that source_rationale is found in the manuscript and that
    group_label/feature_label align within the source_rationale.

    \b
    Examples:
        # Verify and overwrite in place
        mrkr verify markers.json -m manuscript.txt

        # Verify and save to new file
        mrkr verify markers.json -m manuscript.txt -o markers_verified.json
    """
    try:
        manuscript_text = manuscript.read_text(encoding="utf-8")
        records = json.loads(input.read_text(encoding="utf-8"))

        if verbose:
            print(f"\n🔍 Verifying {len(records)} records against {manuscript.name}...")

        records = verify_extractions(manuscript_text, records, verbose=verbose)

        # Count results (only text records get verified)
        text_records = [r for r in records if r.get("source_type") == "text"]
        verified = sum(1 for r in text_records if r.get("_verification", {}).get("all_verified", False))

        out_path = output or input
        out_path.write_text(json.dumps(records, indent=2), encoding="utf-8")

        click.echo(f"\n✅ {verified}/{len(text_records)} text/image extractions fully verified")
        click.echo(f"💾 Saved to: {out_path}")

    except Exception as e:
        click.echo(f"\n❌ Error during verification: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


@cli.command()
@click.argument("input", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--mode",
    required=True,
    type=click.Choice(["celltypes-to-genes", "genes-to-celltypes"]),
    help="Generation mode: celltypes-to-genes or genes-to-celltypes",
)
@click.option(
    "--species",
    "-s",
    required=True,
    help="Species Latin name (e.g., homo_sapiens)",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(path_type=Path),
    help="Output JSON file",
)
@click.option(
    "--metrics",
    type=click.Path(path_type=Path),
    help="Output file for LLM metrics (token usage, costs, timing)",
)
@click.option(
    "--verbose", "-v", is_flag=True, help="Verbose output"
)
def generate(input, mode, species, output, metrics, verbose):
    """
    Generate marker genes or predict cell types using LLM training knowledge.

    Takes an extracted.json file and either generates marker genes for each
    cell type (celltypes-to-genes) or predicts cell types from gene groups
    (genes-to-celltypes). No manuscript or paper context is used.

    \b
    Examples:
        # Generate marker genes for each cell type
        mrkr generate extracted.json --mode celltypes-to-genes -s homo_sapiens -o generated.json

        # Predict cell types from gene groups (anonymized)
        mrkr generate extracted.json --mode genes-to-celltypes -s homo_sapiens -o predicted.json
    """
    # Validate configuration
    try:
        config.validate()
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()

    if not verbose:
        click.echo(f"\n🧬 mrkr generate - {mode}")
        click.echo("=" * 50)

    command_str = " ".join(sys.argv)

    try:
        if mode == "celltypes-to-genes":
            results = generate_celltypes_to_genes(
                input_path=input,
                species=species,
                verbose=verbose,
                metrics_path=metrics,
                command=command_str,
            )
        else:
            results = generate_genes_to_celltypes(
                input_path=input,
                species=species,
                verbose=verbose,
                metrics_path=metrics,
                command=command_str,
            )

        # Save results
        output.write_text(json.dumps(results, indent=2), encoding="utf-8")

        if not verbose:
            click.echo(f"\n✅ Generated {len(results)} marker gene associations")
            click.echo(f"💾 Saved to: {output}")
        else:
            click.echo(f"\n💾 Saved to: {output}")

    except Exception as e:
        click.echo(f"\n❌ Error during generation: {e}", err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    cli()
