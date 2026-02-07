"""Utility functions for mrkr."""

from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd


# Column name variants for flexible matching
COLUMN_VARIANTS = {
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


def find_column(df_cols, variants: List[str]) -> Optional[str]:
    """Find the first matching column (case-insensitive)."""
    for variant in variants:
        for col in df_cols:
            if str(col).lower() == variant.lower():
                return col
    return None


def find_required_columns(df: pd.DataFrame) -> Tuple[bool, Dict[str, str]]:
    """
    Find required columns in DataFrame with flexible naming.

    Returns:
        (found, mapping): found is True if cluster, gene, logfc are found.
                         mapping is dict of key -> actual column name.
    """
    mapping = {}
    for key, variants in COLUMN_VARIANTS.items():
        col = find_column(df.columns, variants)
        if col is not None:
            mapping[key] = col

    # Require cluster, gene, and logfc (p_corr is optional)
    required_keys = ["cluster", "gene", "logfc"]
    found = all(k in mapping for k in required_keys)

    return found, mapping


def load_tabular_with_header_detection(
    file_path: Path,
    sheet_name: Optional[str] = None,
    max_skip: int = 10
) -> Tuple[Optional[pd.DataFrame], Optional[Dict[str, str]]]:
    """
    Load tabular file with automatic header detection.

    Tries different skiprows values (0 to max_skip) to find valid headers.
    For Excel files, optionally loads a specific sheet.

    Returns:
        (df, found_cols): DataFrame and column mapping, or (None, None) if not found.
    """
    file_type = file_path.suffix.lower()

    if file_type == ".xlsx":
        for skip in range(max_skip + 1):
            try:
                if sheet_name:
                    df = pd.read_excel(file_path, sheet_name=sheet_name, skiprows=skip)
                else:
                    df = pd.read_excel(file_path, skiprows=skip)

                found, found_cols = find_required_columns(df)
                if found:
                    return df, found_cols
            except Exception:
                continue

    elif file_type == ".csv":
        for skip in range(max_skip + 1):
            try:
                df = pd.read_csv(file_path, sep=",", skiprows=skip)
                found, found_cols = find_required_columns(df)
                if found:
                    return df, found_cols
            except Exception:
                continue

    elif file_type == ".tsv":
        for skip in range(max_skip + 1):
            try:
                df = pd.read_csv(file_path, sep="\t", skiprows=skip)
                found, found_cols = find_required_columns(df)
                if found:
                    return df, found_cols
            except Exception:
                continue

    return None, None


def get_file_type(file_path: Path) -> str:
    """Determine file type from file extension."""
    suffix = file_path.suffix.lower()

    if suffix == ".xlsx":
        return "xlsx"
    elif suffix == ".csv":
        return "csv"
    elif suffix == ".tsv":
        return "tsv"
    elif suffix == ".txt" or suffix == ".md":
        return "text"
    elif suffix in [".png", ".jpg", ".jpeg"]:
        return "image"
    else:
        raise ValueError(f"Unsupported file type: {suffix}")
