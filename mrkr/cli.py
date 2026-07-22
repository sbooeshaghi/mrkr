"""Command-line interface for source-grounded marker extraction."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import click

from . import __version__
from .claims import (
    ClaimValidationError,
    assert_valid_document,
    validate_document,
)
from .config import config
from .extract import extract_claims
from .ground import ground_document


def write_json_atomic(output: Path, data: object) -> None:
    """Write JSON without exposing a partial output file."""

    output.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=f".{output.name}.", dir=output.parent)
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            json.dump(data, stream, indent=2)
            stream.write("\n")
        os.replace(temporary, output)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def read_json(path: Path) -> object:
    """Read JSON and report a concise command-line error."""

    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise click.ClickException(str(error)) from error


def read_document(path: Path) -> dict:
    """Read one canonical mrkr document."""

    data = read_json(path)
    if not isinstance(data, dict) or not isinstance(data.get("claims"), list):
        raise click.ClickException("expected a canonical mrkr claim or onto document")
    return data


@click.group()
@click.version_option(version=__version__, prog_name="mrkr")
def cli() -> None:
    """Mine source-grounded cell type marker evidence from papers."""


@cli.command()
@click.option(
    "--manuscript",
    "-m",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Source manuscript as Markdown or plain text.",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(path_type=Path),
    help="Output mrkr.claims.v1 JSON file.",
)
@click.option(
    "--source-id",
    help="Stable paper identifier stored in the artifact.",
)
@click.option(
    "--organism",
    required=True,
    help="Organism whose marker claims should be extracted, for example homo_sapiens.",
)
@click.option(
    "--metrics",
    type=click.Path(path_type=Path),
    help="Optional JSON file for model, token, and timing metrics.",
)
@click.option(
    "--response",
    type=click.Path(path_type=Path),
    help="Optional raw Anthropic response JSON retained for audit and parse recovery.",
)
@click.option("--verbose", "-v", is_flag=True, help="Show extraction progress.")
def extract(
    manuscript: Path,
    output: Path,
    source_id: str | None,
    organism: str,
    metrics: Path | None,
    response: Path | None,
    verbose: bool,
) -> None:
    """Extract exact, normalized marker evidence from one manuscript."""

    try:
        config.validate()
        document = extract_claims(
            manuscript_path=manuscript,
            organism=organism,
            source_id=source_id,
            verbose=verbose,
            metrics_path=metrics,
            response_path=response,
            validate=False,
            command="mrkr extract",
        )
        manuscript_text = manuscript.read_text(encoding="utf-8")
        assert_valid_document(document, manuscript_text)
        write_json_atomic(output, document)
    except ClaimValidationError as error:
        stem = output.name.removesuffix(".json")
        write_json_atomic(output.with_name(f"{stem}.rejected.json"), document)
        write_json_atomic(output.with_name(f"{stem}.validation.json"), error.report)
        raise click.ClickException(str(error)) from error
    except click.ClickException:
        raise
    except Exception as error:
        raise click.ClickException(str(error)) from error

    click.echo(f"Extracted {len(document['claims'])} marker claims: {output}")


@cli.command()
@click.argument("input", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(path_type=Path),
    help="Output mrkr.onto.v1 JSON file.",
)
@click.option(
    "--manuscript",
    "-m",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Source manuscript used to verify hashes and exact offsets.",
)
@click.option(
    "--gene-map",
    type=click.Path(exists=True, path_type=Path),
    help="Optional two-column gene symbol to Ensembl map.",
)
@click.option(
    "--organism",
    required=True,
    help="Organism name for gene grounding, for example homo_sapiens.",
)
def ground(
    input: Path,
    output: Path,
    manuscript: Path,
    gene_map: Path | None,
    organism: str,
) -> None:
    """Assign Ensembl, Cell Ontology, and UBERON identifiers."""

    manuscript_text = manuscript.read_text(encoding="utf-8")
    try:
        document = read_document(input)
        grounded = ground_document(
            document,
            organism=organism,
            manuscript_text=manuscript_text,
            gene_map_file=str(gene_map) if gene_map else None,
        )
        write_json_atomic(output, grounded)
    except click.ClickException:
        raise
    except Exception as error:
        raise click.ClickException(str(error)) from error

    terms = [term for claim in grounded["claims"] for term in claim["terms"]]
    resolved = sum(term.get("ontology_term") is not None for term in terms)
    click.echo(
        f"Grounded {len(grounded['claims'])} claims; "
        f"resolved {resolved}/{len(terms)} terms: {output}"
    )


@cli.command("validate")
@click.argument("input", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--manuscript",
    "-m",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Source manuscript used to verify hashes and exact offsets.",
)
@click.option(
    "--report",
    type=click.Path(path_type=Path),
    help="Optional JSON validation report.",
)
def validate_command(input: Path, manuscript: Path, report: Path | None) -> None:
    """Validate a claim or onto document against its source manuscript."""

    manuscript_text = manuscript.read_text(encoding="utf-8")
    try:
        document = read_document(input)
    except (OSError, TypeError, ValueError) as error:
        raise click.ClickException(str(error)) from error

    result = validate_document(document, manuscript_text)
    if report:
        write_json_atomic(report, result)
    for warning in result["warnings"][:10]:
        click.echo(
            f"WARNING {warning['code']} at {warning['path']}: {warning['message']}",
            err=True,
        )
    if result["n_warnings"] > 10:
        click.echo(f"... {result['n_warnings'] - 10} more warnings", err=True)
    if result["valid"]:
        click.echo(
            f"Valid: {result['n_claims']} claims; "
            f"{result['n_warnings']} warnings"
        )
        return
    for error in result["errors"][:10]:
        click.echo(
            f"ERROR {error['code']} at {error['path']}: {error['message']}",
            err=True,
        )
    if result["n_errors"] > 10:
        click.echo(f"... {result['n_errors'] - 10} more errors", err=True)
    raise click.exceptions.Exit(1)


if __name__ == "__main__":
    cli()
