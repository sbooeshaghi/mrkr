"""DEG file processing for marker gene extraction."""

from pathlib import Path
from typing import List, Optional

import pandas as pd

from .utils import load_tabular_with_header_detection, get_file_type


def extract_from_deg_files(
    deg_paths: List[Path],
    species: Optional[str],
    verbose: bool = False,
) -> List[dict]:
    """
    Extract marker genes from DEG table files.

    Handles:
    - Multiple files
    - Excel multi-sheet files
    - CSV/TSV files
    - Flexible column matching
    - Ranking by logFC within cell types

    Args:
        deg_paths: List of DEG file paths
        species: Species name (required for DEG-only extraction)
        verbose: Whether to print progress

    Returns:
        List of evidence dictionaries with DEG metrics
    """
    if not species:
        raise ValueError(
            "Species must be provided when extracting from DEG files. "
            "Use --species option (e.g., --species homo_sapiens)"
        )

    results = []

    for deg_path in deg_paths:
        file_type = get_file_type(deg_path)

        if file_type == "xlsx":
            # Process each sheet as separate data_id
            xl_file = pd.ExcelFile(deg_path)

            if verbose:
                print(f"   📊 Processing {deg_path.name} ({len(xl_file.sheet_names)} sheets)")

            for sheet_name in xl_file.sheet_names:
                try:
                    sheet_results = process_deg_sheet(
                        deg_path, sheet_name, species, verbose
                    )
                    results.extend(sheet_results)
                except ValueError as e:
                    # Skip sheets that don't have valid DEG format
                    if verbose:
                        print(f"      ⚠️  Skipping sheet '{sheet_name}': {str(e).split(chr(10))[0]}")
                    continue

        else:
            # Process CSV/TSV
            if verbose:
                print(f"   📊 Processing {deg_path.name}")

            sheet_results = process_deg_sheet(deg_path, None, species, verbose)
            results.extend(sheet_results)

    return results


def process_deg_sheet(
    file_path: Path,
    sheet_name: Optional[str],
    species: str,
    verbose: bool,
) -> List[dict]:
    """
    Process a single DEG table (file or sheet).

    Args:
        file_path: Path to DEG file
        sheet_name: Sheet name for Excel files (None for CSV/TSV)
        species: Species name
        verbose: Whether to print progress

    Returns:
        List of evidence dictionaries
    """
    # Generate data_id
    file_id = file_path.stem
    if sheet_name:
        data_id = f"{file_id}#{sheet_name}"
    else:
        data_id = file_id

    # Load with header detection
    df, found_cols = load_tabular_with_header_detection(
        file_path, sheet_name=sheet_name
    )

    if df is None:
        # Try to read the file to show what columns are actually there
        try:
            if sheet_name:
                sample_df = pd.read_excel(file_path, sheet_name=sheet_name, nrows=5)
            else:
                sample_df = pd.read_excel(file_path, nrows=5)
            actual_cols = list(sample_df.columns)
            raise ValueError(
                f"Could not find required columns (cluster, gene, logfc) in {data_id}.\n"
                f"Found columns: {actual_cols}\n"
                f"Expected column names (case-insensitive):\n"
                f"  - cluster: cluster, celltype, cell_type, cell-type, group, orig.ident\n"
                f"  - gene: gene, gene_name, genename, gene name\n"
                f"  - logfc: logfc, log_fc, log2fc, avg_log2fc, avg_logFC, etc."
            )
        except Exception:
            raise ValueError(
                f"Could not find required columns (cluster, gene, logfc) in {data_id}. "
                "Please check your DEG table format."
            )

    # Extract column names
    cluster_col = found_cols["cluster"]
    gene_col = found_cols["gene"]
    logfc_col = found_cols.get("logfc")
    pcorr_col = found_cols.get("p_corr")

    # Build evidence records
    records = []
    for _, row in df.iterrows():
        group_label = str(row[cluster_col])
        feature_label = str(row[gene_col])

        # Skip invalid rows
        if pd.isna(group_label) or pd.isna(feature_label):
            continue
        if group_label == "nan" or feature_label == "nan":
            continue

        # Extract metrics
        logfc = None
        if logfc_col and not pd.isna(row[logfc_col]):
            try:
                logfc = float(row[logfc_col])
            except (ValueError, TypeError):
                logfc = None

        pcorr = None
        if pcorr_col and not pd.isna(row[pcorr_col]):
            try:
                pcorr = float(row[pcorr_col])
            except (ValueError, TypeError):
                pcorr = None

        records.append({
            "organism": species,
            "group_label": group_label,
            "group_name": group_label.upper(),
            "group_id": None,
            "feature_label": feature_label,
            "feature_name": feature_label.upper(),
            "feature_id": None,
            "source_type": "deg",
            "source_rationale": "unfiltered",
            "source_id": data_id,
            "data_id": data_id,
            "metrics_pcorr": pcorr,
            "metrics_logfc": logfc,
            "metrics_rank": None,  # Will compute after
        })

    # Compute ranks within each cell type (by logfc descending)
    df_records = pd.DataFrame(records)

    if "metrics_logfc" in df_records.columns and not df_records["metrics_logfc"].isna().all():
        # Rank within each group_name, handling NaNs
        df_records["metrics_rank"] = (
            df_records.groupby("group_name")["metrics_logfc"]
            .rank(ascending=False, method="first", na_option="bottom")
            .astype("Int64")  # Nullable integer type
        )
    else:
        df_records["metrics_rank"] = None

    if verbose:
        n_records = len(df_records)
        n_cell_types = df_records["group_name"].nunique()
        print(f"      ✓ {data_id}: {n_records} markers, {n_cell_types} cell types")

    return df_records.to_dict("records")
