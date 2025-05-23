# tools/heatmap_tool.py

from server import mcp
import scanpy as sc
from utils.plot_utils import save_fig

@mcp.tool()
def plot_rank_genes_heatmap(
    adata_path: str,
    groupby: str = "leiden",
    n_genes: int = 10,
    figdir: str = "figures"
) -> str:
    """
    Plot a heatmap of top-ranked marker genes.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        groupby: Grouping key for marker gene comparison.
        n_genes: Number of top genes per group to plot.
        figdir: Directory to save the heatmap.

    Returns:
        Path to the input AnnData (unchanged, for pipeline compatibility).
    """
    adata = sc.read_h5ad(adata_path)
    fig = sc.pl.rank_genes_groups_heatmap(
        adata, groupby=groupby, n_genes=n_genes, show=False, return_fig=True
    )
    save_fig(fig, f"heatmap_{groupby}.png", figdir)
    return adata_path
