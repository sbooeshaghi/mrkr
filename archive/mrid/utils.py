"""Utility functions for mrid tool"""

import pandas as pd
from typing import Dict, List, Optional, Tuple


def find_column(df_cols, variants: List[str]) -> Optional[str]:
    """Find the first matching column in df for the given list of variants (case-insensitive)."""
    for variant in variants:
        for col in df_cols:
            if str(col).lower() == variant.lower():
                return col
    return None


def find_required_columns(df: pd.DataFrame) -> Tuple[bool, Dict[str, str]]:
    """Find required columns in DataFrame with flexible naming."""
    column_variants = {
        "cluster": [
            "cluster",
            "celltype",
            "cell_type",
            "cell-type",
            "group",
            "orig.ident",
        ],
        "gene": ["gene", "gene_name", "genename", "gene name"],
        "logfc": [
            "logfc",
            "log_fc",
            "log2fc",
            "avg_log2fc",
            "avg_logFC",
            "avg_logfc",
            "log fold change",
            "logfoldchanges",
            "log1p_FC",
        ],
        "p_corr": [
            "p_corr",
            "p_val_adj",
            "pval_adj",
            "pvals_adj",
            "padj",
            "adj_pval",
            "adjusted p-value",
            "pvalue_adj",
            "wilcox.bonferroni",
        ],
    }

    mapping = {}
    for key, variants in column_variants.items():
        col = find_column(df.columns, variants)
        if col is not None:
            mapping[key] = col

    # Only require cluster, gene, and logfc columns
    required_keys = ["cluster", "gene", "logfc"]
    found = all(k in mapping for k in required_keys)

    return found, mapping


def load_excel(
    file_path: str, sheet_name: Optional[str] = None, skiprows: int = 0
) -> pd.DataFrame:
    """Load Excel file with optional sheet name."""
    if sheet_name:
        return pd.read_excel(file_path, sheet_name=sheet_name, skiprows=skiprows)
    else:
        # Load the first readable sheet
        excel_file = pd.ExcelFile(file_path)

        for sheet in excel_file.sheet_names:
            try:
                df = pd.read_excel(file_path, sheet_name=sheet, skiprows=skiprows)
                if not df.empty:
                    return df
            except Exception:
                # Skip sheets that can't be read
                continue

        raise ValueError(f"No readable sheets found in {file_path}")


def load_csv_tsv(file_path: str, sep: str = ",", skiprows: int = 0) -> pd.DataFrame:
    """Load CSV or TSV file."""
    return pd.read_csv(file_path, sep=sep, skiprows=skiprows)


def load_tabular(file_path: str, file_type: str, skiprows: int = 0) -> pd.DataFrame:
    """Load a tabular file (csv, tsv, xlsx) into a pandas DataFrame."""
    if file_type == "xlsx":
        return load_excel(file_path, skiprows=skiprows)
    elif file_type == "csv":
        return load_csv_tsv(file_path, sep=",", skiprows=skiprows)
    elif file_type == "tsv":
        return load_csv_tsv(file_path, sep="\t", skiprows=skiprows)
    else:
        raise ValueError(f"Unsupported file type: {file_type}")


def load_tabular_with_header_detection(
    file_path: str, file_type: str, max_skip: int = 10
) -> Tuple[Optional[pd.DataFrame], Optional[Dict[str, str]]]:
    """Try loading the tabular file with skiprows=0..max_skip, returning the first DataFrame where all required columns are found."""

    if file_type == "xlsx":
        # For Excel files, try each sheet individually
        try:
            excel_file = pd.ExcelFile(file_path)
            for sheet_name in excel_file.sheet_names:
                for skip in range(0, max_skip + 1):
                    try:
                        df = pd.read_excel(
                            file_path, sheet_name=sheet_name, skiprows=skip
                        )
                        found, found_cols = find_required_columns(df)
                        if found:
                            return df, found_cols
                    except Exception:
                        continue
        except Exception:
            pass
    else:
        # For other file types, use the original approach
        for skip in range(0, max_skip + 1):
            try:
                df = load_tabular(file_path, file_type, skiprows=skip)
                found, found_cols = find_required_columns(df)
                if found:
                    return df, found_cols
            except Exception:
                continue

    return None, None


def get_file_id(file_name: str) -> str:
    """Extract file_id from file_name (handles Excel sheet references)."""
    if "#" in file_name:
        return file_name
    return file_name


def get_file_type(file_path: str) -> str:
    """Determine file type from file path."""
    if file_path.endswith(".xlsx"):
        return "xlsx"
    elif file_path.endswith(".csv"):
        return "csv"
    elif file_path.endswith(".tsv"):
        return "tsv"
    elif file_path.endswith(".txt"):
        return "txt"
    elif file_path.lower().endswith((".png", ".jpg", ".jpeg")):
        return "image"
    else:
        raise ValueError(f"Unsupported file type: {file_path}")
