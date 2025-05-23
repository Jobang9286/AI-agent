# tools/velocity_tool.py

from server import mcp
import scanpy as sc
import scvelo as scv
from utils.plot_utils import save_fig
import os

@mcp.tool()
def compute_rna_velocity(
    adata_path: str,
    loom_file: str = None,
    mode: str = "stochastic",
    figdir: str = "figures"
) -> str:
    """
    Perform RNA velocity analysis using scVelo.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        loom_file: Path to loom file with spliced and unspliced counts.
                  If None, will try to use existing data in adata.
        mode: Velocity mode ('deterministic', 'stochastic', or 'dynamical').
        figdir: Directory to save velocity visualizations.

    Returns:
        Path to the .h5ad file with velocity results.
    """
    # Load the data
    adata = sc.read_h5ad(adata_path)
    
    # Check if we need to load velocity data from loom file
    if loom_file and os.path.exists(loom_file):
        # Load loom file and merge with existing annotations
        ldata = scv.read(loom_file, cache=True)
        adata = scv.utils.merge(adata, ldata)
    
    # Process the data for velocity analysis
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    
    # Compute velocity
    if mode == "deterministic":
        scv.tl.velocity(adata, mode='deterministic')
    elif mode == "dynamical":
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode='dynamical')
    else:  # Default to stochastic
        scv.tl.velocity(adata, mode='stochastic')
    
    # Calculate velocity graph and pseudotime
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_pseudotime(adata)
    
    # Save plots
    # Velocity stream plot
    scv.pl.velocity_embedding_stream(adata, basis='umap', color='velocity_pseudotime', show=False, save=False)
    save_fig(None, "velocity_stream.png", figdir, plt_figure=True)
    
    # Velocity pseudotime plot
    scv.pl.scatter(adata, color='velocity_pseudotime', show=False, save=False)
    save_fig(None, "velocity_pseudotime.png", figdir, plt_figure=True)
    
    # Save results
    out_path = adata_path.replace(".h5ad", "_velocity.h5ad")
    adata.write(out_path)
    return out_path