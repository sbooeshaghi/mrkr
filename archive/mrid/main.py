"""Main CLI for mrid tool"""

import argparse
import os
import sys

from .commands import (
    extract_deg_command,
    # extract_image_command,
    extract_text_command,
    init_command,
    map_command,
)
from .models import InferenceModel


def read_context_file(context_path: str) -> str:
    """Read context from file if path is provided, otherwise return as-is."""
    if not context_path:
        return ""

    # Check if context looks like a file path
    if os.path.exists(context_path) and os.path.isfile(context_path):
        try:
            with open(context_path, "r", encoding="utf-8") as f:
                return f.read().strip()
        except Exception as e:
            print(f"Warning: Could not read context file {context_path}: {e}")
            return context_path
    else:
        # Not a file path, return as-is
        return context_path


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="mrid - Cell-type marker gene extraction tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  mrid init deg_file.xlsx -o spec.tsv
  mrid extract -f deg deg_file.xlsx -s homo_sapiens -x spec.tsv -o evidence.json
  mrid extract -f text manuscript.txt paper.txt -s homo_sapiens -x spec.tsv -o evidence.json
  mrid extract -f image figure1.png figure2.png -s homo_sapiens -x spec.tsv -o evidence.json
  mrid map evidence.json -o mapped.json
  mrid map evidence.json -o mapped.json --gene-map custom_mapping.tsv
        """,
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Init command
    init_parser = subparsers.add_parser(
        "init", help="Create standardized spec from DEG files"
    )
    init_parser.add_argument(
        "input_files", nargs="+", help="Input DEG files (Excel, CSV, TSV)"
    )
    init_parser.add_argument("-o", "--output", required=True, help="Output spec file")
    init_parser.add_argument(
        "--context", default="", help="Study context or path to context file"
    )
    init_parser.add_argument(
        "--model",
        choices=[m.name for m in InferenceModel],
        default=InferenceModel.OPENAI_GPT4o.name,
        help="LLM model to use",
    )

    # Extract command
    extract_parser = subparsers.add_parser(
        "extract", help="Extract marker genes from files"
    )
    extract_parser.add_argument(
        "-f",
        "--format",
        choices=["deg", "text", "image"],
        required=True,
        help="Input format",
    )
    extract_parser.add_argument("input_files", nargs="+", help="Input files")
    extract_parser.add_argument("-s", "--species", required=True, help="Species name")
    extract_parser.add_argument("-x", "--spec", required=True, help="Spec file")
    extract_parser.add_argument("-o", "--output", required=True, help="Output file")
    extract_parser.add_argument(
        "--context", default="", help="Study context or path to context file"
    )
    extract_parser.add_argument(
        "--model",
        choices=[m.name for m in InferenceModel],
        default=InferenceModel.OPENAI_GPT4o.name,
        help="LLM model to use",
    )

    # Map command
    map_parser = subparsers.add_parser("map", help="Map gene names to standardized IDs")
    map_parser.add_argument("input_file", help="Input evidence file")
    map_parser.add_argument("-o", "--output", required=True, help="Output file")
    map_parser.add_argument(
        "--gene-map",
        help="Gene mapping file (default: uses package's mappable_gnames.tsv)",
    )
    map_parser.add_argument(
        "--context", default="", help="Study context or path to context file"
    )
    map_parser.add_argument(
        "--model",
        choices=[m.name for m in InferenceModel],
        default=InferenceModel.OPENAI_GPT4o.name,
        help="LLM model to use",
    )

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    # Parse model
    try:
        model = InferenceModel[args.model]
    except KeyError:
        print(f"Error: Invalid model '{args.model}'")
        sys.exit(1)

    # Read context from file if it's a file path
    context = read_context_file(args.context)

    try:
        if args.command == "init":
            init_command(args.input_files, args.output, model, context)

        elif args.command == "extract":
            if args.format == "deg":
                extract_deg_command(
                    args.input_files,
                    args.species,
                    args.spec,
                    args.output,
                    context,
                    model,
                )
            elif args.format == "text":
                extract_text_command(
                    args.input_files,
                    args.species,
                    args.spec,
                    args.output,
                    context,
                    model,
                )
            # elif args.format == "image":
            #     extract_image_command(
            #         args.input_files,
            #         args.species,
            #         args.spec,
            #         args.output,
            #         context,
            #         model,
            #     )

        elif args.command == "map":
            map_command(args.input_file, args.output, args.gene_map, context, model)

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
