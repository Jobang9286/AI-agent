# tools/umap_tool.py

from server import mcp
import scanpy as sc
from utils.plot_utils import save_fig

@mcp.tool()
def compute_umap(adata_path: str, figdir: str = "figures", color: str = "leiden") -> str:
    """
    Compute UMAP projection and save a UMAP plot.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        figdir: Directory to save the UMAP plot.
        color: Column in adata.obs to color the UMAP by.

    Returns:
        Path to the .h5ad file with UMAP coordinates computed.
    """
    adata = sc.read_h5ad(adata_path)
    sc.tl.umap(adata)

    fig = sc.pl.umap(adata, color=color, show=False, return_fig=True)
    save_fig(fig, f"umap_{color}.png", figdir)

    out_path = adata_path.replace(".h5ad", "_umap.h5ad")
    adata.write(out_path)
    return out_path
