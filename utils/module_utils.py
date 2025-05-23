# utils/module_utils.py

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Optional, Union, Tuple
import os
from utils.plot_utils import save_fig

def identify_gene_modules_nmf(adata: ad.AnnData, n_components: int = 20, 
                             use_hvg: bool = True, random_state: int = 42):
    """
    Identify gene modules using Non-negative Matrix Factorization (NMF).
    
    Args:
        adata: AnnData object
        n_components: Number of modules to identify
        use_hvg: Whether to use highly variable genes
        random_state: Random seed for reproducibility
        
    Returns:
        Tuple of (cell embeddings, gene loadings)
    """
    from sklearn.decomposition import NMF
    
    # Use HVG if requested and available
    if use_hvg and 'highly_variable' in adata.var:
        X = adata[:, adata.var.highly_variable].X
    else:
        X = adata.X
    
    # Convert to dense if sparse
    if hasattr(X, 'toarray'):
        X = X.toarray()
    
    # Run NMF
    model = NMF(n_components=n_components, random_state=random_state)
    cell_embeddings = model.fit_transform(X)  # Cells x Modules
    gene_loadings = model.components_  # Modules x Genes
    
    return cell_embeddings, gene_loadings
def store_gene_modules(adata: ad.AnnData, cell_embeddings: np.ndarray, 
                      gene_loadings: np.ndarray, gene_names: List[str], 
                      n_top_genes: int = 100, score_cells: bool = True):
    """
    Store gene modules in AnnData object.
    
    Args:
        adata: AnnData object
        cell_embeddings: Cell embedding matrix (cells x modules)
        gene_loadings: Gene loading matrix (modules x genes)
        gene_names: List of gene names corresponding to columns in gene_loadings
        n_top_genes: Number of top genes to include per module
        score_cells: Whether to compute module scores for cells
        
    Returns:
        Updated AnnData object
    """
    # Store cell embeddings
    adata.obsm['module_embeddings'] = cell_embeddings
    
    # Create dataframe of gene loadings
    gene_loadings_df = pd.DataFrame(gene_loadings, columns=gene_names)
    
    # Store in uns
    adata.uns['gene_modules'] = {}
    adata.uns['gene_modules']['loadings'] = gene_loadings_df
    adata.uns['gene_modules']['n_modules'] = cell_embeddings.shape[1]
    
    # Store top genes per module and compute scores
    module_genes = {}
    for i in range(cell_embeddings.shape[1]):
        module_name = f'module_{i+1}'
        
        # Get top genes for this module
        top_genes = gene_loadings_df.iloc[i].nlargest(n_top_genes).index.tolist()
        module_genes[module_name] = top_genes
        
        # Score cells by module activity
        if score_cells:
            sc.tl.score_genes(adata, gene_list=top_genes, score_name=module_name)
    
    adata.uns['gene_modules']['top_genes'] = module_genes
    return adata
def plot_module_heatmap(adata: ad.AnnData, n_modules: int = 6, n_genes_per_module: int = 10,
                       groupby: str = 'leiden', figdir: str = "figures", 
                       filename: str = "module_heatmap.png"):
    """
    Plot heatmap of module genes.
    
    Args:
        adata: AnnData object with gene modules
        n_modules: Number of modules to include in the plot
        n_genes_per_module: Number of top genes to include per module
        groupby: Column in adata.obs to group cells by
        figdir: Directory to save figure
        filename: Name of output file
    """
    if 'gene_modules' not in adata.uns or 'top_genes' not in adata.uns['gene_modules']:
        print("No gene modules found in adata.uns")
        return
    
    # Get top genes for each module
    all_genes = []
    for i in range(min(n_modules, adata.uns['gene_modules']['n_modules'])):
        module_name = f'module_{i+1}'
        if module_name in adata.uns['gene_modules']['top_genes']:
            module_genes = adata.uns['gene_modules']['top_genes'][module_name]
            all_genes.extend(module_genes[:n_genes_per_module])
    
    if not all_genes:
        print("No genes found in modules")
        return
    
    # Filter to unique genes
    unique_genes = list(dict.fromkeys(all_genes))
    
    # Create AnnData object for heatmap
    adata_heatmap = adata[:, unique_genes].copy()
    
    # Add module annotation to variables
    adata_heatmap.var['module'] = 'Other'
    for i in range(min(n_modules, adata.uns['gene_modules']['n_modules'])):
        module_name = f'module_{i+1}'
        if module_name in adata.uns['gene_modules']['top_genes']:
            module_genes = adata.uns['gene_modules']['top_genes'][module_name][:n_genes_per_module]
            for gene in module_genes:
                if gene in adata_heatmap.var_names:
                    adata_heatmap.var.at[gene, 'module'] = module_name
    
    # Plot heatmap
    sc.pl.heatmap(adata_heatmap, var_names=unique_genes, groupby=groupby, 
                 var_group_positions=[(0, len(unique_genes))],
                 var_group_labels=['Module genes'],
                 show=False, save=False)
    
    save_fig(None, filename, figdir, plt_figure=True)
def plot_module_scores(adata: ad.AnnData, basis: str = 'umap', 
                      n_modules: int = 6, ncols: int = 3,
                      figdir: str = "figures", 
                      filename: str = "module_scores.png"):
    """
    Plot module scores on embedding.
    
    Args:
        adata: AnnData object with module scores
        basis: Basis for plotting (e.g., 'umap', 'tsne')
        n_modules: Number of modules to plot
        ncols: Number of columns in the plot grid
        figdir: Directory to save figure
        filename: Name of output file
    """
    if f'X_{basis}' not in adata.obsm:
        print(f"Warning: {basis} embedding not found in adata.obsm")
        return
    
    # Check which module scores are available
    module_columns = []
    for i in range(n_modules):
        module_name = f'module_{i+1}'
        if module_name in adata.obs:
            module_columns.append(module_name)
    
    if not module_columns:
        print("No module scores found in adata.obs")
        return
    
    # Plot module scores
    fig = sc.pl.embedding(adata, basis=basis, color=module_columns, 
                         ncols=ncols, show=False, return_fig=True)
    save_fig(fig, filename, figdir)
def module_gene_network(adata: ad.AnnData, module_id: int, 
                       n_top_genes: int = 30, 
                       corr_threshold: float = 0.3,
                       figdir: str = "figures", 
                       filename: str = None):
    """
    Create a gene correlation network for a specific module.
    
    Args:
        adata: AnnData object with gene modules
        module_id: Module ID (1-based)
        n_top_genes: Number of top genes to include
        corr_threshold: Correlation threshold for drawing edges
        figdir: Directory to save figure
        filename: Name of output file (defaults to module_{module_id}_network.png)
    """
    try:
        import networkx as nx
    except ImportError:
        print("networkx package not installed. Please install it with 'pip install networkx'")
        return
    
    if 'gene_modules' not in adata.uns or 'top_genes' not in adata.uns['gene_modules']:
        print("No gene modules found in adata.uns")
        return
    
    module_name = f'module_{module_id}'
    if module_name not in adata.uns['gene_modules']['top_genes']:
        print(f"Module {module_id} not found")
        return
    
    # Get top genes for this module
    top_genes = adata.uns['gene_modules']['top_genes'][module_name][:n_top_genes]
    
    # Calculate gene correlation matrix
    adata_sub = adata[:, top_genes].copy()
    if hasattr(adata_sub.X, 'toarray'):
        expr_matrix = pd.DataFrame(adata_sub.X.toarray(), columns=top_genes)
    else:
        expr_matrix = pd.DataFrame(adata_sub.X, columns=top_genes)
    
    corr_matrix = expr_matrix.corr()    
    # Create network
    G = nx.Graph()
    
    # Add nodes
    for gene in top_genes:
        G.add_node(gene)
    
    # Add edges based on correlation threshold
    for i, gene1 in enumerate(top_genes):
        for gene2 in top_genes[i+1:]:
            corr = corr_matrix.loc[gene1, gene2]
            if abs(corr) >= corr_threshold:
                G.add_edge(gene1, gene2, weight=abs(corr))
    
    # Draw network
    plt.figure(figsize=(12, 12))
    
    # Calculate node size based on gene importance in the module
    gene_importance = adata.uns['gene_modules']['loadings'].loc[module_id-1, top_genes]
    max_importance = gene_importance.max()
    min_importance = gene_importance.min()
    
    node_size = {gene: 100 + 500 * (gene_importance[gene] - min_importance) / (max_importance - min_importance) 
                for gene in top_genes}
    
    # Use spring layout
    pos = nx.spring_layout(G, k=0.3, iterations=50)
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, 
                          node_size=[node_size[gene] for gene in G.nodes],
                          node_color='skyblue')
    
    # Draw edges with width proportional to correlation
    for (gene1, gene2, data) in G.edges(data=True):
        width = data['weight'] * 2
        nx.draw_networkx_edges(G, pos, edgelist=[(gene1, gene2)], width=width, alpha=0.7)
    
    # Draw labels
    nx.draw_networkx_labels(G, pos, font_size=10)
    
    plt.title(f"Gene Network for {module_name}")
    plt.axis('off')
    plt.tight_layout()
    
    # Save figure
    if filename is None:
        filename = f"module_{module_id}_network.png"
    
    os.makedirs(figdir, exist_ok=True)
    plt.savefig(f"{figdir}/{filename}", dpi=300, bbox_inches='tight')
    plt.close()
    
    return G
def identify_gene_modules_pca(adata: ad.AnnData, n_components: int = 20,
                             use_hvg: bool = True, random_state: int = 42):
    """
    Identify gene modules using Principal Component Analysis (PCA).
    
    Args:
        adata: AnnData object
        n_components: Number of components to identify
        use_hvg: Whether to use highly variable genes
        random_state: Random seed for reproducibility
        
    Returns:
        Tuple of (cell embeddings, gene loadings)
    """
    from sklearn.decomposition import PCA
    
    # Use HVG if requested and available
    if use_hvg and 'highly_variable' in adata.var:
        X = adata[:, adata.var.highly_variable].X
    else:
        X = adata.X
    
    # Convert to dense if sparse
    if hasattr(X, 'toarray'):
        X = X.toarray()
    
    # Run PCA
    pca = PCA(n_components=n_components, random_state=random_state)
    cell_embeddings = pca.fit_transform(X)  # Cells x Components
    gene_loadings = pca.components_  # Components x Genes
    
    return cell_embeddings, gene_loadings
def compare_modules(adata: ad.AnnData, groupby: str = 'leiden',
                   figdir: str = "figures", 
                   filename: str = "module_vs_clusters.png"):
    """
    Compare identified modules with cell clusters.
    
    Args:
        adata: AnnData object with gene modules and clusters
        groupby: Column in adata.obs indicating clusters
        figdir: Directory to save figure
        filename: Name of output file
    """
    if 'gene_modules' not in adata.uns:
        print("No gene modules found in adata.uns")
        return
    
    # Get module scores
    module_columns = []
    for i in range(adata.uns['gene_modules']['n_modules']):
        module_name = f'module_{i+1}'
        if module_name in adata.obs:
            module_columns.append(module_name)
    
    if not module_columns:
        print("No module scores found in adata.obs")
        return
    
    # Calculate average module score per cluster
    module_scores_by_cluster = pd.DataFrame(index=adata.obs[groupby].unique())
    
    for module in module_columns:
        module_scores_by_cluster[module] = adata.obs.groupby(groupby)[module].mean()
    
    # Plot heatmap
    plt.figure(figsize=(12, 8))
    im = plt.imshow(module_scores_by_cluster, cmap='viridis')
    plt.colorbar(im, label='Mean Module Score')
    
    # Add labels
    plt.yticks(range(len(module_scores_by_cluster.index)), module_scores_by_cluster.index)
    plt.xticks(range(len(module_scores_by_cluster.columns)), module_scores_by_cluster.columns, 
              rotation=45, ha='right')
    
    plt.title('Module Scores Across Cell Clusters')
    plt.tight_layout()
    
    # Save figure
    os.makedirs(figdir, exist_ok=True)
    plt.savefig(f"{figdir}/{filename}", dpi=300, bbox_inches='tight')
    plt.close()
    
    return module_scores_by_cluster