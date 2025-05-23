# tools/pca_tool.py

from server import mcp
import scanpy as sc

@mcp.tool()
def compute_pca(adata_path: str, n_comps: int = 50) -> str:
    """
    Compute PCA and store the result in the AnnData object.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        n_comps: Number of principal components to compute.

    Returns:
        Path to the .h5ad file with PCA results.
    """
    adata = sc.read_h5ad(adata_path)
    sc.tl.pca(adata, n_comps=n_comps)

    out_path = adata_path.replace(".h5ad", "_pca.h5ad")
    adata.write(out_path)
    return out_path
