import pandas as pd
import scanpy as sc
import scipy.sparse


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
