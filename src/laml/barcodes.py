"""Helpers for TCGA barcodes.

Barcode example: TCGA-AB-2934-03A-01R-0000-00
  - characters 1-12  -> patient ID   (TCGA-AB-2934)
  - characters 14-15 -> sample type  (01-09 tumor, 10-19 normal)
"""

from __future__ import annotations

import pandas as pd


def patient_id(barcode: str) -> str:
    """Return the 12-character patient ID from a TCGA barcode."""
    return barcode[:12]


def sample_type_code(barcode: str) -> str:
    """Return the two-digit sample type code (characters 14-15)."""
    return barcode[13:15]


def is_tumor(barcode: str) -> bool:
    """True for tumor sample types (01-09), False for normals (10-19)."""
    code = int(sample_type_code(barcode))
    return 1 <= code <= 9


def keep_tumor_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only tumor-sample columns of a features x samples matrix."""
    return df.loc[:, [c for c in df.columns if is_tumor(c)]]


def one_sample_per_patient(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse to one column per patient (first occurrence, sorted by barcode).

    Replicates are rare but possible; the rule is explicit so results are
    reproducible. Columns are renamed to the 12-character patient ID.
    """
    cols = sorted(df.columns)
    seen: dict[str, str] = {}
    for c in cols:
        seen.setdefault(patient_id(c), c)
    out = df.loc[:, list(seen.values())].copy()
    out.columns = list(seen.keys())
    return out
