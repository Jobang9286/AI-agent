# tools/pseudotime_tool.py

from server import mcp
import scanpy as sc
from utils.plot_utils import save_fig

@mcp.tool()
def compute_pseudotime(
    adata_path: str,
    root_cells: str = None, 
    n_dcs: int = 15,
    use_diffmap: bool = True,
    figdir: str = "figures"
) -> str:
    """
    Perform pseudotime analysis using diffusion pseudotime.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        root_cells: Comma-separated list of cell IDs to use as root cells.
                   If None, the root will be identified automatically.
        n_dcs: Number of diffusion components to use.
        use_diffmap: Whether to compute diffusion map first.
        figdir: Directory to save pseudotime visualizations.

    Returns:
        Path to the .h5ad file with pseudotime results.
    """
    adata = sc.read_h5ad(adata_path)
    
    # Check if neighbors have been computed
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata)
    
    # Compute diffusion map if requested
    if use_diffmap:
        sc.tl.diffmap(adata, n_comps=n_dcs)
    
    # Set up root cells
    root = None
    if root_cells:
        root_cell_ids = [cell_id.strip() for cell_id in root_cells.split(',')]
        if all(cell_id in adata.obs_names for cell_id in root_cell_ids):
            root = adata.obs_names.get_indexer(root_cell_ids)[0]
    
    # Compute diffusion pseudotime
    sc.tl.dpt(adata, n_dcs=n_dcs, root=root)
    
    # Save plots
    if 'X_umap' in adata.obsm:
        # Plot pseudotime on UMAP
        fig = sc.pl.umap(adata, color='dpt_pseudotime', show=False, return_fig=True)
        save_fig(fig, "pseudotime_umap.png", figdir)
    
    if 'X_diffmap' in adata.obsm:
        # Plot pseudotime on diffusion map
        fig = sc.pl.diffmap(adata, color='dpt_pseudotime', show=False, return_fig=True)
        save_fig(fig, "pseudotime_diffmap.png", figdir)
    
    # Save results
    out_path = adata_path.replace(".h5ad", "_pseudotime.h5ad")
    adata.write(out_path)
    return out_path