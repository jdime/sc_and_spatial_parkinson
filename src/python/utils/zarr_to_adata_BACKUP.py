import os
import sys
import anndata as ad
import scanpy as sc

from anndata.experimental import read_dispatched, write_dispatched, read_elem
import dask.array as da
import zarr
from collections.abc import Mapping
from typing import Any

import warnings

warnings.filterwarnings('ignore')

def read_dask(store):
    f = zarr.open(store, mode="r")

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