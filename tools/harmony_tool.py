# tools/harmony_tool.py
from server import mcp
import scanpy as sc
import harmonypy
from utils.plot_utils import save_fig
import matplotlib.pyplot as plt
from pathlib import Path

@mcp.tool()
def compute_harmony(adata_path: str, batch_key: str, n_pcs: int = 30, 
                   theta: float = 2.0, figdir: str = "figures") -> str:
    """
    Perform Harmony integration on PCA embedding for batch correction.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        batch_key: Key in adata.obs for batch annotation.
        n_pcs: Number of principal components to use.
        theta: Diversity clustering penalty parameter.
        figdir: Directory to save the harmony embedding plot.

    Returns:
        Path to the .h5ad file with harmony-corrected embedding.
    """
    adata = sc.read_h5ad(adata_path)
    
    # Ensure PCA has been computed
    if 'X_pca' not in adata.obsm:
        sc.pp.pca(adata, n_comps=n_pcs)
    
    # Run Harmony
    pca_matrix = adata.obsm['X_pca'][:, :n_pcs]
    batch_labels = adata.obs[batch_key].values
    
    harmony_model = harmonypy.RunHarmony(
        pca_matrix,
        batch_labels,
        theta=theta
    )
    
    # Store the corrected embedding
    adata.obsm['X_harmony'] = harmony_model.Z_corr
    
    # Create UMAP with harmony embeddings
    sc.pp.neighbors(adata, use_rep='X_harmony')
    sc.tl.umap(adata)
    
    # Plot and save
    fig = sc.pl.umap(adata, color=batch_key, show=False, return_fig=True)
    save_fig(fig, f"harmony_umap_{batch_key}.png", figdir)
    
    out_path = adata_path.replace(".h5ad", "_harmony.h5ad")
    adata.write(out_path)
    return out_path