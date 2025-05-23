# tools/clustering_tool.py

from server import mcp
import scanpy as sc

@mcp.tool()
def leiden_clustering(adata_path: str, resolution: float = 1.0) -> str:
    """
    Perform Leiden clustering on the neighborhood graph.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        resolution: Resolution parameter for the Leiden algorithm.

    Returns:
        Path to the .h5ad file with clustering labels.
    """
    adata = sc.read_h5ad(adata_path)
    sc.tl.leiden(adata, resolution=resolution)

    out_path = adata_path.replace(".h5ad", "_leiden.h5ad")
    adata.write(out_path)
    return out_path
