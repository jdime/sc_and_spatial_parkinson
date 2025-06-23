import scanpy as sc
import pandas as pd
import logging

from typing import List, Optional

logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)


def add_metadata_to_adata(
    adata: sc.AnnData,
    metadata_df: pd.DataFrame,
    join: List[str],
    columns_to_add: List[str],
) -> sc.AnnData:
    """Adds metadata to AnnData object `obs` slot

    Args:
        adata (sc.AnnData): AnnData object with `obs` slot to add metadata
        metadata_df (pd.DataFrame): Dataframe containing barcodes in rows and metadata in columns
        join (List[str]): Indicates key(s) to join the AnnData `obs` and metadata (e.g. `barcode` or `study_name + sample_name`)
        columns_to_add (List[str]): Columns to add to AnnData `obs`. These columns must be present in `metadata_df`

    Raises:
        ValueError: Mismatch in shape of original vs. merged AnnData `obs`
        ValueError: Merge did not preserve order of barcodes (i.e. the left join messed up the order of barcodes)

    Returns:
        sc.AnnData: AnnData object with cell metadata added from dataframe
    """
    original_adata_columns = list(adata.obs.columns)
    logger.info(f"Adding metadata to adata, joining on {join}, keeping columns {columns_to_add}...")
    new_columns = list(set(columns_to_add).intersection(set(metadata_df.columns)))

    for col in new_columns:
        if col in adata.obs.columns and col not in join:
            logger.info(f"{col} in adata, dropping")
            adata.obs.drop(columns=col, inplace=True, axis=1)

    merged_obs = adata.obs.merge(metadata_df, how="left", on=join)[original_adata_columns + new_columns]
    # merged_obs = merged_obs.iloc[:, :-4] ### Javier's patch because it's duplicating these columns

    merged_obs["index"] = merged_obs.barcode
    merged_obs.set_index("index", inplace=True)

    if merged_obs.shape[0] != adata.obs.shape[0]:
        raise ValueError("Mismatch in shape of original adata obs and merged obs")
    if not (merged_obs.barcode.values == adata.obs.barcode.values).all():
        raise ValueError("The merge did not preserve order of barcodes")

    adata.obs = merged_obs

    return adata
