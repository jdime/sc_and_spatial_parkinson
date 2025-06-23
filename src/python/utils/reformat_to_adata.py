import os
import sys
import anndata as ad
import scanpy as sc
import pandas as pd
import scipy.sparse

### tested with zarr==2.18.7 and dask==2025.4.0
from anndata.experimental import read_dispatched, write_dispatched, read_elem
import dask.array as da 
import zarr
from collections.abc import Mapping
from typing import Any

import warnings
warnings.filterwarnings('ignore')

def convert_zarr_to_adata(
    zarr_dir
) -> sc.AnnData:
    """Creates an AnnData object out of a zarr directory

    Args:
        zarr_dir (str): Path to a zarr file with subfolders: layers, obs, obsm, obsp, uns, var, varm, varp, X

    Returns:
        sc.AnnData: AnnData object
    """
    
    f = zarr.open(zarr_dir, mode="r")

    def callback(func, elem_name: str, elem, iospec):
        if iospec.encoding_type in (
            "dataframe",
            "csr_matrix",
            "csc_matrix",
            "awkward-array",
        ):
            # Preventing recursing inside of these types
            return read_elem(elem)
        elif iospec.encoding_type == "array":
            return da.from_zarr(elem)
        else:
            return func(elem)
    adata = read_dispatched(f, callback=callback)

    return adata


def convert_mtx_to_adata(
    matrix: scipy.sparse,
    barcodes_df: pd.DataFrame,
    features_df: pd.DataFrame,
    barcodes_prefix: str,
) -> sc.AnnData:
    """Creates an AnnData object out of a matrix, barcodes_df and features_df. Optionally add a prefix to barcodes

    Args:
        matrix (scipy.sparse): Sparse array corresponding to matrix
        barcodes_df (pd.DataFrame): Dataframe corresponding to barcodes
        features_df (pd.DataFrame): Dataframe corresponding to features
        barcodes_prefix: str | None (default: None)

    Returns:
        sc.AnnData: AnnData object
    """
    ### Make unique gene names (equivalent to scanpy.read_10x_mtx(.., make_unique=True)
    if len(features_df["gene"].unique()) != len(features_df["gene"]):
        gb = features_df.groupby([features_df.gene])
        features_df['gene'] = features_df.gene + [str(-a) if a else '' for a in gb.cumcount()]

    ### Add prefix to barcodes
    if barcodes_prefix != None:
        barcodes_df["barcode"] = barcodes_prefix + barcodes_df["barcode"]

    assert len(barcodes_df["barcode"].unique()) == len(barcodes_df["barcode"]), "Barcodes are not unique."
    assert len(features_df["gene"].unique()) == len(features_df["gene"]), "Genes are not unique."
    adata = sc.AnnData(X=matrix.T.tocsr(), obs=barcodes_df, var=features_df)
    adata.obs["index"] = adata.obs.barcode
    adata.obs.set_index("index", inplace=True)

    return adata

    