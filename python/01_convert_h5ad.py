#!/usr/bin/env python3
"""
01_convert_h5ad.py

GSE243292 Alzheimer's disease snRNA-seq
H5AD -> Seurat-compatible sparse Matrix Market conversion.

Python is used ONLY for data conversion/preparation.
All biological analysis is performed later in R using Seurat.

Default input:
D:/nsg consultaion/data/raw/GSE243292/GSE243292_ADsnRNAseq_GEO.h5ad

Default output:
D:/nsg consultaion/data/processed

Outputs:
    matrix.mtx.gz       sparse genes x cells expression matrix
    features.tsv.gz     gene ID and gene name
    barcodes.tsv.gz     cell IDs
    cell_metadata.tsv   AnnData obs metadata
    gene_metadata.tsv   AnnData var metadata
    conversion_report.json

Memory safety:
    The script refuses dense matrices. This avoids creating a large dense
    copy that could exhaust a Windows laptop with ~12 GB RAM.
"""

from __future__ import annotations

import argparse
import gzip
import json
import platform
import shutil
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.io import mmwrite


DEFAULT_INPUT = Path(
    r"D:\nsg consultaion\data\raw\GSE243292\GSE243292_ADsnRNAseq_GEO.h5ad"
)

DEFAULT_OUTPUT = Path(
    r"D:\nsg consultaion\data\processed"
)


def parse_args() -> argparse.Namespace:
    """Read optional input/output paths and matrix selection."""
    parser = argparse.ArgumentParser(
        description="Convert GSE243292 H5AD to sparse Seurat-compatible files."
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Input .h5ad file."
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output directory."
    )

    parser.add_argument(
        "--matrix-source",
        choices=("auto", "counts", "x", "raw"),
        default="auto",
        help=(
            "Expression matrix source. 'auto' uses layers['counts'] "
            "if present, otherwise X. 'raw' uses raw.X."
        )
    )

    return parser.parse_args()


def clean_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Convert object columns to strings for robust TSV export."""
    result = df.copy()

    for column in result.columns:
        if pd.api.types.is_object_dtype(result[column]):
            result[column] = result[column].astype(str)

    return result


def select_matrix(adata: ad.AnnData, source: str):
    """
    Select the expression matrix.

    auto:
        layers['counts'] if available, otherwise X
    counts:
        layers['counts']
    x:
        X
    raw:
        raw.X
    """

    if source == "counts":
        if "counts" not in adata.layers:
            raise KeyError("layers['counts'] is not present in this H5AD.")
        return adata.layers["counts"], "layers/counts"

    if source == "x":
        return adata.X, "X"

    if source == "raw":
        if adata.raw is None:
            raise ValueError("adata.raw is not available in this H5AD.")
        return adata.raw.X, "raw.X"

    if "counts" in adata.layers:
        return adata.layers["counts"], "layers/counts"

    return adata.X, "X"


def write_mtx_gzip(matrix, output_file: Path) -> None:
    """Write a sparse Matrix Market file and gzip it."""
    temp_file = output_file.with_suffix("")

    mmwrite(str(temp_file), matrix)

    with temp_file.open("rb") as source:
        with gzip.open(output_file, "wb", compresslevel=6) as destination:
            shutil.copyfileobj(
                source,
                destination,
                length=1024 * 1024
            )

    temp_file.unlink(missing_ok=True)


def main() -> int:
    """Run the complete conversion."""
    args = parse_args()

    input_path = args.input.expanduser().resolve()
    output_dir = args.output.expanduser().resolve()

    print("=" * 70)
    print("GSE243292 H5AD -> Seurat conversion")
    print("=" * 70)

    # ------------------------------------------------------------
    # 1. Validate input
    # ------------------------------------------------------------
    print("\n[1/8] Checking input file...")

    if not input_path.exists():
        raise FileNotFoundError(
            f"\nH5AD file not found:\n{input_path}\n"
        )

    if input_path.suffix.lower() != ".h5ad":
        raise ValueError("Input must be an .h5ad file.")

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Input : {input_path}")
    print(f"Output: {output_dir}")

    # ------------------------------------------------------------
    # 2. Read AnnData
    # ------------------------------------------------------------
    print("\n[2/8] Reading H5AD...")

    # AnnData officially supports reading H5AD files with read_h5ad().
    adata = ad.read_h5ad(input_path)

    print(adata)
    print(f"Cells: {adata.n_obs:,}")
    print(f"Genes: {adata.n_vars:,}")

    # ------------------------------------------------------------
    # 3. Select expression matrix
    # ------------------------------------------------------------
    print("\n[3/8] Selecting expression matrix...")

    matrix, matrix_source = select_matrix(
        adata,
        args.matrix_source
    )

    print(f"Matrix source: {matrix_source}")

    # Do not allow an accidental dense conversion.
    if not sparse.issparse(matrix):
        raise MemoryError(
            "\nThe selected matrix is DENSE.\n"
            "The script stopped intentionally for memory safety.\n\n"
            "For a 12 GB RAM laptop, inspect the H5AD and try:\n"
            "  --matrix-source counts\n"
            "or\n"
            "  --matrix-source x\n"
        )

    print("Matrix storage: sparse")
    print(f"Non-zero entries: {matrix.nnz:,}")

    # ------------------------------------------------------------
    # 4. Orient as genes x cells
    # ------------------------------------------------------------
    print("\n[4/8] Preparing genes x cells matrix...")

    # AnnData stores observations/cells x variables/genes.
    # Seurat ReadMtx expects features/genes x cells.
    matrix_gene_cell = matrix.tocsc().T.tocsc()

    print(
        f"Export dimensions: "
        f"{matrix_gene_cell.shape[0]:,} genes x "
        f"{matrix_gene_cell.shape[1]:,} cells"
    )

    # ------------------------------------------------------------
    # 5. Export gene and cell identifiers/metadata
    # ------------------------------------------------------------
    print("\n[5/8] Writing feature and cell metadata...")

    gene_ids = pd.Index(
        adata.var_names.astype(str),
        name="gene_id"
    )

    var = adata.var.copy()

    gene_name_column = None

    for candidate in (
        "gene_name",
        "gene_symbol",
        "symbol",
        "gene"
    ):
        if candidate in var.columns:
            gene_name_column = candidate
            break

    if gene_name_column is not None:
        gene_names = (
            var[gene_name_column]
            .astype(str)
            .reset_index(drop=True)
        )
    else:
        gene_names = pd.Series(
            gene_ids.astype(str),
            name="gene_name"
        )

    # Seurat feature file:
    # column 1 = gene ID
    # column 2 = gene name
    features = pd.DataFrame({
        "gene_id": gene_ids.astype(str),
        "gene_name": gene_names.astype(str)
    })

    features.to_csv(
        output_dir / "features.tsv.gz",
        sep="\t",
        header=False,
        index=False,
        compression="gzip"
    )

    # Preserve complete gene metadata.
    gene_metadata = clean_metadata(
        var.reset_index(drop=True)
    )

    gene_metadata.insert(
        0,
        "gene_id",
        gene_ids.astype(str)
    )

    if "gene_name" not in gene_metadata.columns:
        gene_metadata.insert(
            1,
            "gene_name",
            gene_names.astype(str)
        )

    gene_metadata.to_csv(
        output_dir / "gene_metadata.tsv",
        sep="\t",
        index=False
    )

    # Cell IDs / barcodes.
    cell_ids = pd.Index(
        adata.obs_names.astype(str),
        name="cell_id"
    )

    pd.DataFrame({
        "cell_id": cell_ids.astype(str)
    }).to_csv(
        output_dir / "barcodes.tsv.gz",
        sep="\t",
        header=False,
        index=False,
        compression="gzip"
    )

    # Complete AnnData cell metadata.
    cell_metadata = clean_metadata(
        adata.obs.copy()
    )

    cell_metadata.index = cell_metadata.index.astype(str)

    cell_metadata.insert(
        0,
        "cell_id",
        cell_metadata.index
    )

    cell_metadata.to_csv(
        output_dir / "cell_metadata.tsv",
        sep="\t",
        index=False
    )

    # ------------------------------------------------------------
    # 6. Write sparse Matrix Market
    # ------------------------------------------------------------
    print("\n[6/8] Writing sparse matrix.mtx.gz...")

    write_mtx_gzip(
        matrix_gene_cell,
        output_dir / "matrix.mtx.gz"
    )

    # ------------------------------------------------------------
    # 7. Write conversion report
    # ------------------------------------------------------------
    print("\n[7/8] Writing conversion report...")

    report = {
        "project": "GSE243292 Alzheimer's disease snRNA-seq",
        "input_file": str(input_path),
        "output_directory": str(output_dir),
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "matrix_source": matrix_source,
        "exported_rows_genes": int(matrix_gene_cell.shape[0]),
        "exported_columns_cells": int(matrix_gene_cell.shape[1]),
        "sparse": True,
        "nonzero_entries": int(matrix_gene_cell.nnz),
        "gene_name_source": (
            gene_name_column
            if gene_name_column is not None
            else "var_names"
        ),
        "python_version": sys.version,
        "platform": platform.platform(),
        "anndata_version": ad.__version__,
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "scipy_version": scipy.__version__,
        "analysis_note": (
            "Python is used only for conversion/preprocessing. "
            "Main biological analysis is performed in R using Seurat."
        )
    }

    with (
        output_dir / "conversion_report.json"
    ).open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)

    # ------------------------------------------------------------
    # 8. Final summary
    # ------------------------------------------------------------
    print("\n[8/8] Conversion completed.")
    print("=" * 70)
    print("Generated files:")
    print("  matrix.mtx.gz")
    print("  features.tsv.gz")
    print("  barcodes.tsv.gz")
    print("  cell_metadata.tsv")
    print("  gene_metadata.tsv")
    print("  conversion_report.json")
    print("=" * 70)

    print(
        "\nNext step in RStudio:\n"
        'source("R/01_load_converted_data.R")'
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
