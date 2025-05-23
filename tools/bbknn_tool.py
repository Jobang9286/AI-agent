# tools/bbknn_tool.py
from server import mcp
import scanpy as sc
import bbknn
import os

@mcp.tool()
def compute_bbknn(adata_path: str, batch_key: str, neighbors_within_batch: int = 3, 
                 n_pcs: int = 50, trim: float = 0) -> str:
    """
    Compute batch balanced k-nearest neighbors (BBKNN) for scRNA-seq integration.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        batch_key: Key in adata.obs for batch annotation.
        neighbors_within_batch: Number of neighbors to consider within each batch.
        n_pcs: Number of principal components to use.
        trim: Trim outliers with this fraction of batches.

    Returns:
        Path to the .h5ad file with the batch-corrected neighborhood graph.
    """
    adata = sc.read_h5ad(adata_path)
    
    bbknn.bbknn(
        adata, 
        batch_key=batch_key,
        neighbors_within_batch=neighbors_within_batch,
        n_pcs=n_pcs,
        trim=trim,
        copy=False
    )
    
    out_path = adata_path.replace(".h5ad", "_bbknn.h5ad")
    adata.write(out_path)
    return out_path