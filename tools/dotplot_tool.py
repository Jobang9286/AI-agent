# tools/dotplot_tool.py

from server import mcp
import scanpy as sc
from utils.plot_utils import save_fig

@mcp.tool()
def plot_dotplot(
    adata_path: str,
    marker_genes: list[str],
    groupby: str = "leiden",
    figdir: str = "figures"
) -> str:
    """
    Generate a dotplot for the given marker genes across clusters.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        marker_genes: List of marker gene names.
        groupby: Column to group cells by.
        figdir: Directory to save the figure.

    Returns:
        Path to the input AnnData (unchanged, for pipeline compatibility).
    """
    adata = sc.read_h5ad(adata_path)
    fig = sc.pl.dotplot(adata, marker_genes, groupby=groupby, show=False, return_fig=True)
    save_fig(fig, f"dotplot_{groupby}.png", figdir)
    return adata_path
