"""Command implementations for mrid tool"""

import json
import os
from pathlib import Path
from typing import List, Optional

import pandas as pd

from .llm import (
    extract_from_text_with_combined_enums,
    infer_cell_names,

)
from .models import Evidence, InferenceModel, SpecRecord
from .utils import (
    get_file_id,
    get_file_type,
    load_tabular_with_header_detection,
    find_required_columns,
)


def init_command(
    input_files: List[str],
    output_file: str,
    model: InferenceModel = InferenceModel.OPENAI_GPT4o,
    context: str = "",
) -> None:
    """Create standardized spec from DEG files."""
    spec_records = []
    all_unique_labels = set()  # Collect all unique labels across all files

    for file_path in input_files:
        if not os.path.exists(file_path):
            print(f"Warning: File {file_path} does not exist, skipping")
            continue

        file_name = os.path.basename(file_path)
        file_type = get_file_type(file_path)

        if file_type == "xlsx":
            try:
                excel_file = pd.ExcelFile(file_path)
                for sheet_name in excel_file.sheet_names:
                    df, found_cols = None, None
                    for skip in range(0, 11):  # Try up to 10 skip rows
                        try:
                            df = pd.read_excel(
                                file_path, sheet_name=sheet_name, skiprows=skip
                            )
                            found, found_cols = find_required_columns(df)
                            if found:
                                break
                        except Exception:
                            continue

                    if df is None or found_cols is None:
                        print(f"Sheet '{sheet_name}' does not have required columns")
                        continue

                    required_keys = ["cluster", "gene", "logfc"]
                    missing_keys = [
                        key for key in required_keys if key not in found_cols
                    ]
                    if missing_keys:
                        print(
                            f"Sheet '{sheet_name}' does not have columns: {', '.join(missing_keys)}"
                        )
                        continue

                    # Extract unique cell labels from this sheet
                    cell_labels = df[found_cols["cluster"]].unique().tolist()
                    all_unique_labels.update(cell_labels)  # Add to global set

                    # Create spec records for this sheet (without LLM inference yet)
                    for label in cell_labels:
                        spec_records.append(
                            SpecRecord(
                                file_id=f"{os.path.splitext(file_name)[0]}#{sheet_name}",
                                file_name=file_name,
                                file_type=file_type,
                                file_uri=os.path.abspath(file_path),
                                data_id=f"{os.path.splitext(file_name)[0]}#{sheet_name}",
                                data_type="deg",
                                group_label=str(label),
                                group_name="",  # Will be filled after LLM inference
                                group_id="",
                            )
                        )

            except Exception as e:
                print(f"Error processing Excel file {file_path}: {e}")
                continue

        else:
            # Handle CSV/TSV files
            df, found_cols = load_tabular_with_header_detection(file_path, file_type)
            if df is None or found_cols is None:
                print(f"File '{file_name}' does not have required columns")
                continue

            required_keys = ["cluster", "gene", "logfc"]
            missing_keys = [key for key in required_keys if key not in found_cols]
            if missing_keys:
                print(
                    f"File '{file_name}' does not have columns: {', '.join(missing_keys)}"
                )
                continue

            # Extract unique cell labels from this file
            cell_labels = df[found_cols["cluster"]].unique().tolist()
            all_unique_labels.update(cell_labels)  # Add to global set

            # Create spec records for this file (without LLM inference yet)
            for label in cell_labels:
                spec_records.append(
                    SpecRecord(
                        file_id=get_file_id(file_name),
                        file_name=file_name,
                        file_type=file_type,
                        file_uri=os.path.abspath(file_path),
                        data_id=get_file_id(file_name),
                        data_type="deg",
                        group_label=str(label),
                        group_name="",  # Will be filled after LLM inference
                        group_id="",
                    )
                )

    # Now do batch LLM inference for all unique labels
    print(f"Found {len(all_unique_labels)} unique labels across all files")
    if all_unique_labels:
        try:
            standardized = infer_cell_names(list(all_unique_labels), model, context)

            # Create mapping from original label to standardized name
            label_to_standardized = {}
            for item in standardized:
                label_to_standardized[item["original_label"]] = item[
                    "standardized_name"
                ]

            # Update all spec records with standardized names
            for record in spec_records:
                record.group_name = label_to_standardized.get(
                    record.group_label, record.group_label
                )

        except Exception as e:
            print(f"Warning: Error during LLM inference: {e}")
            # If LLM fails, use original labels
            for record in spec_records:
                record.group_name = record.group_label

    # Write spec file
    spec_df = pd.DataFrame([record.model_dump() for record in spec_records])
    spec_df.to_csv(output_file, sep="\t", index=False)
    print(f"Created spec file: {output_file}")
    print(f"Spec file shape: {spec_df.shape}")
    print(f"Spec file columns: {list(spec_df.columns)}")


def extract_deg_command(
    deg_files: List[str],
    species: str,
    spec_file: str,
    output_file: str,
    context: str = "",
    model: InferenceModel = InferenceModel.OPENAI_GPT4o,
) -> None:
    """Extract marker genes from DEG tables."""
    # Load spec file
    spec_df = pd.read_csv(spec_file, sep="\t")
    evidence_list = []

    for file_path in deg_files:
        file_type = get_file_type(file_path)
        file_name = Path(file_path).name

        if file_type not in ["xlsx", "csv", "tsv"]:
            print(f"Skipping {file_path} - unsupported file type")
            continue

        # Get all file IDs that match this file name
        # For Excel files, this will include sheet-specific IDs (filename#sheetname)
        # For other files, this will be just the filename
        file_specs = spec_df[spec_df["file_name"] == file_name]

        if file_specs.empty:
            print(f"No spec found for file: {file_name}")
            continue

        # Load and find columns
        df, found_cols = load_tabular_with_header_detection(file_path, file_type)
        if df is None or found_cols is None:
            print(f"Skipping {file_path} - could not find required columns")
            continue

        # Process each unique file_id for this file
        unique_file_ids = file_specs["file_id"].unique()

        for file_id in unique_file_ids:
            print(f"Processing file_id: {file_id}")

            # Get spec records for this file_id
            file_id_specs = file_specs[file_specs["file_id"] == file_id]

            # Load data for this file_id
            if file_type == "xlsx" and "#" in file_id:
                # Extract sheet name from file_id (format: filename#sheetname)
                sheet_name = file_id.split("#", 1)[1]
                try:
                    # Try to load this specific sheet
                    sheet_df, sheet_found_cols = None, None
                    for skip in range(0, 11):  # Try up to 10 skip rows
                        try:
                            sheet_df = pd.read_excel(
                                file_path, sheet_name=sheet_name, skiprows=skip
                            )
                            found, sheet_found_cols = find_required_columns(sheet_df)
                            if found:
                                break
                        except Exception:
                            continue

                    if sheet_df is None or sheet_found_cols is None:
                        print(f"Sheet '{sheet_name}' does not have required columns")
                        continue

                    df = sheet_df
                    found_cols = sheet_found_cols

                except Exception as e:
                    print(f"Error loading sheet '{sheet_name}': {e}")
                    continue
            else:
                # Use the already loaded data for non-Excel files
                df = df
                found_cols = found_cols

            # Merge data with spec records using pandas operations
            df["cluster"] = df[found_cols["cluster"]].astype(str)
            df["gene"] = df[found_cols["gene"]].astype(str)
            df["logfc"] = (
                df[found_cols["logfc"]].astype(float)
                if found_cols.get("logfc")
                else 0.0
            )
            df["p_corr"] = (
                df[found_cols["p_corr"]].astype(float)
                if found_cols.get("p_corr")
                else -1
            )

            # Create a mapping from cluster to spec info
            cluster_to_spec = file_id_specs.set_index("group_label")[
                ["group_name", "group_id"]
            ].to_dict("index")

            # Filter data to only include clusters that exist in spec
            valid_clusters = set(cluster_to_spec.keys())
            df_filtered = df[df["cluster"].isin(valid_clusters)]

            # Create evidence records using vectorized operations
            for _, row in df_filtered.iterrows():
                cluster_val = row["cluster"]
                spec_info = cluster_to_spec[cluster_val]

                evidence = Evidence(
                    organism=species,
                    group_label=cluster_val,
                    group_name=str(spec_info["group_name"]),
                    group_id=str(spec_info["group_id"])
                    if pd.notna(spec_info["group_id"])
                    else "",
                    feature_label=row["gene"],
                    source_id=file_id,
                    data_id=file_id,
                    metrics_pcorr=row["p_corr"],
                    metrics_logfc=row["logfc"],
                )
                evidence_list.append(evidence)

    # Calculate ranks
    evidence_df = pd.DataFrame([e.model_dump() for e in evidence_list])
    if not evidence_df.empty:
        evidence_df["metrics_rank"] = evidence_df.groupby("group_name")[
            "metrics_logfc"
        ].rank(ascending=False)
        evidence_list = [
            Evidence.model_validate(row.to_dict()) for _, row in evidence_df.iterrows()
        ]

    # Write JSON output
    with open(output_file, "w") as f:
        json.dump([e.model_dump() for e in evidence_list], f, indent=2)

    print(f"Extracted evidence to: {output_file}")


def extract_text_command(
    text_files: List[str],
    species: str,
    spec_file: str,
    output_file: str,
    context: str = "",
    model: InferenceModel = InferenceModel.OPENAI_GPT4o,
) -> None:
    """Extract marker genes from text files using enum-based cell type mapping."""
    # Load spec file
    spec_df = pd.read_csv(spec_file, sep="\t")
    evidence_list = []

    for file_path in text_files:
        if not file_path.endswith(".txt"):
            print(f"Skipping {file_path} - not a text file")
            continue

        print(f"Processing text file: {file_path}")

        # Read text content
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                text_content = f.read()
        except Exception as e:
            print(f"Error reading file {file_path}: {e}")
            continue

        # Extract associations using LLM with combined enum-based mapping (single pass)
        extractions = extract_from_text_with_combined_enums(
            text_content, spec_df, context, model
        )

        if not extractions:
            print(f"No extractions found in {file_path}")
            continue

        # Filter out invalid group_name and data_id pairs
        # valid_extractions = []
        # for extraction in extractions:
        #     group_name = extraction.get("group_name", "")
        #     data_id = extraction.get("data_id", "")

        #     # Check if this combination exists in the spec
        #     matching_spec = spec_df[
        #         (spec_df["group_name"] == group_name) & (spec_df["data_id"] == data_id)
        #     ]

        #     if not matching_spec.empty:
        #         valid_extractions.append(extraction)
        #     else:
        #         print(
        #             f"Filtering out invalid pair: group_name='{group_name}', data_id='{data_id}'"
        #         )

        # extractions = valid_extractions

        # Create evidence records
        print(len(extractions))
        print(extractions)
        for extraction in extractions:
            group_label = extraction.get("group_label", "")
            feature_label = extraction.get("feature_label", "")
            group_name = extraction.get("group_name", "")
            data_id = extraction.get("data_id", "")
            source_rationale = extraction.get("source_rationale", "")

            if not group_label or not feature_label or not group_name or not data_id:
                continue

            # Find matching spec record to get group_id (data_id is already from extraction)
            matching_spec = spec_df[(spec_df["group_name"] == group_name.value)]
            if matching_spec.empty:
                print(
                    f"No matching spec record found for group_name: {group_name}, data_id: {data_id}"
                )
                continue

            spec_row = matching_spec.iloc[0]

            evidence = Evidence(
                organism=species,
                group_label=group_label,
                group_name=group_name,
                group_id=str(spec_row["group_id"])
                if pd.notna(spec_row["group_id"])
                else "",
                feature_label=feature_label,
                feature_name=feature_label,
                feature_id=None,
                source_type="text",
                source_rationale=source_rationale,
                source_id=Path(file_path).name,
                data_id=data_id,
                metrics_pcorr=None,
                metrics_logfc=None,
                metrics_rank=None,
            )
            evidence_list.append(evidence)

    # Convert to DataFrame for easier manipulation
    if evidence_list:
        evidence_df = pd.DataFrame([e.model_dump() for e in evidence_list])
        print(
            f"Extracted {len(evidence_df)} associations from {len(text_files)} text files"
        )
    else:
        evidence_df = pd.DataFrame()
        print("No evidence extracted from text files")

    # Write JSON output
    with open(output_file, "w") as f:
        json.dump([e.model_dump() for e in evidence_list], f, indent=2)

    print(f"Extracted evidence to: {output_file}")


# def extract_image_command(
#     image_files: List[str],
#     species: str,
#     spec_file: str,
#     output_file: str,
#     context: str = "",
#     model: InferenceModel = InferenceModel.OPENAI_GPT4o,
# ) -> None:
#     """Extract marker genes from image files."""
#     # Load spec file
#     spec_df = pd.read_csv(spec_file, sep="\t")
#     evidence_list = []

#     for file_path in image_files:
#         if Path(file_path).suffix.lower() not in [".png", ".jpg", ".jpeg"]:
#             print(f"Skipping {file_path} - not an image file")
#             continue

#         print(f"Processing image file: {file_path}")

#         # Extract associations using vision LLM
#         try:
#             extractions = extract_from_image(file_path, context, model)
#         except Exception as e:
#             print(f"Error processing image {file_path}: {e}")
#             continue

#         if not extractions:
#             print(f"No extractions found in {file_path}")
#             continue

#         # Create evidence records
#         for extraction in extractions:
#             cell_type = extraction.get("cell_type", "")
#             gene = extraction.get("gene", "")
#             extraction_context = extraction.get("context", "")
#             data_id = extraction.get("data_id", "")

#             if not cell_type or not gene:
#                 continue

#             # Find matching spec record - try exact match first, then fallback
#             matching_spec = spec_df[
#                 spec_df["group_name"].str.lower() == cell_type.lower()
#             ]
#             if matching_spec.empty:
#                 # Use first spec record as fallback
#                 matching_spec = spec_df.head(1)

#             if not matching_spec.empty:
#                 spec_row = matching_spec.iloc[0]

#                 evidence = Evidence(
#                     organism=species,
#                     group_label=cell_type,
#                     group_name=spec_row["group_name"],
#                     group_id=str(spec_row["group_id"])
#                     if pd.notna(spec_row["group_id"])
#                     else "",
#                     feature_label=gene,
#                     feature_name=None,
#                     feature_id=None,
#                     source_type="image",
#                     source_rationale=extraction_context,
#                     source_id=Path(file_path).name,
#                     data_id=data_id or spec_row["data_id"],
#                     metrics_pcorr=None,
#                     metrics_logfc=None,
#                     metrics_rank=None,
#                 )
#                 evidence_list.append(evidence)

#     # Convert to DataFrame for easier manipulation
#     if evidence_list:
#         evidence_df = pd.DataFrame([e.model_dump() for e in evidence_list])
#         print(
#             f"Extracted {len(evidence_df)} associations from {len(image_files)} image files"
#         )
#     else:
#         evidence_df = pd.DataFrame()
#         print("No evidence extracted from image files")

#     # Write JSON output
#     with open(output_file, "w") as f:
#         json.dump([e.model_dump() for e in evidence_list], f, indent=2)

#     print(f"Extracted evidence to: {output_file}")


def map_command(
    extracted_json: str,
    output_file: str,
    gene_map_file: Optional[str] = None,
    context: str = "",
    model: InferenceModel = InferenceModel.OPENAI_GPT4o,
) -> None:
    """Map gene names to standardized IDs."""
    # Use default mappable_gnames.tsv from package data if no file provided
    if gene_map_file is None:
        package_dir = Path(__file__).parent.parent
        gene_map_file = str(package_dir / "data" / "mappable_gnames.tsv")

    # Load mappable gene names (only where mappable is True)
    df = (
        pd.read_csv(gene_map_file, sep="\t", keep_default_na=False, na_values=[])
        .query("mappable")
        .set_index("gname")
    )

    # Load extracted evidence with pandas
    ext = pd.read_json(extracted_json)

    # Map feature_name to feature_id using the mappable genes
    n = ext.feature_name.map(df["gid"])
    u = ext.feature_name.str.upper().map(df["gid"])
    ext["feature_id"] = n.combine_first(u)

    # For unmapped genes, use LLM fallback if feature_id is still null
    unmapped_mask = ext["feature_id"].isna()
    if unmapped_mask.any():
        unmapped_genes = ext.loc[unmapped_mask, "feature_name"].unique().tolist()
        print(f"Found {len(unmapped_genes)} unmapped genes, using LLM fallback")
        print(unmapped_genes)



    # Write updated JSON
    ext.to_json(output_file, orient="records", indent=2)

    print(f"Mapped evidence to: {output_file}")
