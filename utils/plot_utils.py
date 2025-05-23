# utils/plot_utils.py

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import os
import numpy as np

def plot_umap(adata: ad.AnnData, figdir: str = "figures", color: str = "leiden"):
    fig = sc.pl.umap(adata, color=color, show=False, return_fig=True)
    save_fig(fig, f"umap_{color}.png", figdir)

def plot_umap_cluster(adata: ad.AnnData, figdir: str = "figures", color: str = "leiden"):
    fig = sc.pl.umap(adata, color=color, show=False, return_fig=True)
    save_fig(fig, f"umap_cluster_{color}.png", figdir)

def save_fig(fig, filename: str, figdir: str = "figures", dpi: int = 300, close: bool = True, plt_figure: bool = False):
    """
    Save a matplotlib figure to the specified directory.

    Args:
        fig: Matplotlib figure object. If plt_figure=True, can be None to use current plt figure.
        filename: Name of the file to save (e.g., "umap.png").
        figdir: Directory to save the figure.
        dpi: Image resolution.
        close: Whether to close the figure after saving.
        plt_figure: Whether to use plt.gcf() instead of provided fig.
    """
    os.makedirs(figdir, exist_ok=True)
    path = os.path.join(figdir, filename)
    
    if plt_figure:
        plt.savefig(path, dpi=dpi, bbox_inches='tight')
    else:
        if fig is not None:
            fig.savefig(path, dpi=dpi, bbox_inches='tight')
        else:
            plt.savefig(path, dpi=dpi, bbox_inches='tight')
    
    if close:
        plt.close()

def plot_marker_genes_heatmap(adata: ad.AnnData, marker_genes: list, groupby: str = "leiden", 
                             figdir: str = "figures", n_cells: int = 20):
    """
    Plot heatmap of marker genes across clusters.
    
    Args:
        adata: AnnData object
        marker_genes: List of marker genes to plot
        groupby: Column in adata.obs to group cells by
        figdir: Directory to save figure
        n_cells: Number of cells per group to plot
    """
    sc.pl.heatmap(adata, marker_genes, groupby=groupby, 
                 n_cells=n_cells, show=False, save=False)
    save_fig(None, f"marker_genes_heatmap_{groupby}.png", figdir, plt_figure=True)

def plot_dotplot(adata: ad.AnnData, marker_genes: list, groupby: str = "leiden", 
                figdir: str = "figures"):
    """
    Plot dot plot of marker genes across clusters.
    
    Args:
        adata: AnnData object
        marker_genes: List of marker genes to plot
        groupby: Column in adata.obs to group cells by
        figdir: Directory to save figure
    """
    sc.pl.dotplot(adata, marker_genes, groupby=groupby, 
                 show=False, save=False)
    save_fig(None, f"dotplot_{groupby}.png", figdir, plt_figure=True)

def plot_violin(adata: ad.AnnData, genes: list, groupby: str = "leiden", 
               figdir: str = "figures", rotation: float = 90):
    """
    Plot violin plot of genes across groups.
    
    Args:
        adata: AnnData object
        genes: List of genes to plot
        groupby: Column in adata.obs to group cells by
        figdir: Directory to save figure
        rotation: Rotation angle for x-axis labels
    """
    sc.pl.violin(adata, genes, groupby=groupby, rotation=rotation,
                show=False, save=False)
    save_fig(None, f"violin_{groupby}.png", figdir, plt_figure=True)

def plot_trajectory(adata: ad.AnnData, basis: str = "umap", color: str = "dpt_pseudotime", 
                   figdir: str = "figures"):
    """
    Plot trajectory on a dimensionality reduction embedding.
    
    Args:
        adata: AnnData object with trajectory information
        basis: Basis for plotting (e.g., 'umap', 'tsne', 'diffmap')
        color: Column in adata.obs to color by (e.g., 'dpt_pseudotime')
        figdir: Directory to save figure
    """
    if f'X_{basis}' not in adata.obsm:
        print(f"Warning: {basis} embedding not found in adata.obsm")
        return
    
    fig = sc.pl.embedding(adata, basis=basis, color=color, show=False, return_fig=True)
    save_fig(fig, f"trajectory_{basis}_{color}.png", figdir)

def plot_gene_expression(adata: ad.AnnData, genes: list, basis: str = "umap", 
                        ncols: int = 3, figdir: str = "figures"):
    """
    Plot gene expression on dimensionality reduction embedding.
    
    Args:
        adata: AnnData object
        genes: List of genes to plot
        basis: Basis for plotting (e.g., 'umap', 'tsne')
        ncols: Number of columns in the plot grid
        figdir: Directory to save figure
    """
    if f'X_{basis}' not in adata.obsm:
        print(f"Warning: {basis} embedding not found in adata.obsm")
        return
    
    # Filter genes to those present in the data
    genes_to_plot = [gene for gene in genes if gene in adata.var_names]
    
    if not genes_to_plot:
        print("Warning: None of the specified genes found in the data")
        return
    
    fig = sc.pl.embedding(adata, basis=basis, color=genes_to_plot, ncols=ncols, 
                         show=False, return_fig=True)
    save_fig(fig, f"gene_expression_{basis}.png", figdir)

def plot_barplot(values, labels, title: str, xlabel: str, ylabel: str, 
                figdir: str = "figures", filename: str = "barplot.png", 
                figsize=(10, 6), color: str = 'skyblue', horizontal: bool = False):
    """
    Create a customizable bar plot.
    
    Args:
        values: Array of values for the bars
        labels: Labels for the bars
        title: Plot title
        xlabel: X-axis label
        ylabel: Y-axis label
        figdir: Directory to save figure
        filename: Name of the output file
        figsize: Size of the figure as a tuple (width, height)
        color: Color of the bars
        horizontal: Whether to create a horizontal bar plot
    """
    plt.figure(figsize=figsize)
    
    if horizontal:
        plt.barh(range(len(values)), values, color=color)
        plt.yticks(range(len(labels)), labels)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    else:
        plt.bar(range(len(values)), values, color=color)
        plt.xticks(range(len(labels)), labels, rotation=45, ha='right')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    
    plt.title(title)
    plt.tight_layout()
    
    save_fig(None, filename, figdir, plt_figure=True)
