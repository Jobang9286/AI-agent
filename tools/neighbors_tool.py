# tools/neighbors_tool.py

from server import mcp
import scanpy as sc

@mcp.tool()
def compute_neighbors(adata_path: str, n_neighbors: int = 10, n_pcs: int = 40) -> str:
    """
    Compute the neighborhood graph of cells using PCA representation.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        n_neighbors: Number of neighbors for KNN graph.
        n_pcs: Number of principal components to use.

    Returns:
        Path to the .h5ad file with neighbor graph computed.
    """
    adata = sc.read_h5ad(adata_path)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    out_path = adata_path.replace(".h5ad", "_neighbors.h5ad")
    adata.write(out_path)
    return out_path
