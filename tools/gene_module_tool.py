# tools/gene_module_tool.py

from server import mcp
import scanpy as sc
import pandas as pd
import numpy as np
from utils.plot_utils import save_fig

@mcp.tool()
def analyze_gene_modules(
    adata_path: str,
    n_modules: int = 20,
    method: str = "nmf",
    figdir: str = "figures"
) -> str:
    """
    Identify gene modules and regulatory networks from single-cell data.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        n_modules: Number of gene modules to identify.
        method: Method to identify gene modules ('nmf', 'pca', 'ica').
        figdir: Directory to save gene module visualizations.

    Returns:
        Path to the .h5ad file with gene module results.
    """
    adata = sc.read_h5ad(adata_path)
    
    # Store the full data for later
    adata_full = adata.copy()
    
    # Use highly variable genes if available, otherwise use all genes
    if 'highly_variable' in adata.var:
        adata = adata[:, adata.var.highly_variable]
    
    # Extract expression matrix
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    
    # Identify gene modules using the specified method
    if method == "nmf":
        from sklearn.decomposition import NMF
        model = NMF(n_components=n_modules, random_state=42)
        W = model.fit_transform(X)  # Cell embeddings
        H = model.components_       # Gene embeddings
        
        # Store module scores and gene loadings
        adata_full.obsm['NMF'] = W
        
        # Create a dataframe of gene loadings
        gene_loadings = pd.DataFrame(H, columns=adata.var_names)
        
        # Assign genes to modules
        for i in range(n_modules):
            module_name = f'module_{i+1}'
            # Get top genes for module
            top_genes = gene_loadings.iloc[i].nlargest(100).index.tolist()
            
            # Score cells by module activity
            sc.tl.score_genes(adata_full, gene_list=top_genes, score_name=module_name)
            
            # Store module genes in adata.uns
            if 'gene_modules' not in adata_full.uns:
                adata_full.uns['gene_modules'] = {}
            adata_full.uns['gene_modules'][module_name] = top_genes        
        # Visualize module scores
        if 'X_umap' in adata_full.obsm:
            module_columns = [f'module_{i+1}' for i in range(min(6, n_modules))]
            fig = sc.pl.umap(adata_full, color=module_columns, show=False, return_fig=True)
            save_fig(fig, "gene_modules_umap.png", figdir)
        
        # Generate a heatmap of top genes per module
        top_module_genes = []
        for i in range(min(6, n_modules)):
            module_name = f'module_{i+1}'
            top_genes = gene_loadings.iloc[i].nlargest(10).index.tolist()
            top_module_genes.extend(top_genes)
        
        if len(top_module_genes) > 0:
            adata_heatmap = adata_full[:, top_module_genes]
            # Add module annotation to var
            adata_heatmap.var['module'] = 'Other'
            for i in range(min(6, n_modules)):
                module_name = f'module_{i+1}'
                top_genes = gene_loadings.iloc[i].nlargest(10).index.tolist()
                for gene in top_genes:
                    if gene in adata_heatmap.var_names:
                        adata_heatmap.var.at[gene, 'module'] = module_name
            
            # Plot heatmap
            sc.pl.heatmap(adata_heatmap, var_names=top_module_genes, groupby='module', 
                         dendrogram=True, show=False, save=f"{figdir}/module_heatmap.png")
    
    # Save results
    out_path = adata_path.replace(".h5ad", "_modules.h5ad")
    adata_full.write(out_path)
    return out_path