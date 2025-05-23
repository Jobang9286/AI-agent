# utils/velocity_utils.py

import scvelo as scv
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from typing import Optional, List, Union
import os
from utils.plot_utils import save_fig

def setup_velocity(adata, n_pcs: int = 30, n_neighbors: int = 30, 
                  min_shared_counts: int = 20, n_top_genes: int = 2000):
    """
    Set up the AnnData object for velocity analysis.

    Args:
        adata: AnnData object with spliced and unspliced counts.
        n_pcs: Number of principal components to use.
        n_neighbors: Number of neighbors for graph.
        min_shared_counts: Minimum number of shared counts.
        n_top_genes: Number of top genes to use.

    Returns:
        AnnData object ready for velocity computation.
    """
    scv.pp.filter_and_normalize(adata, min_shared_counts=min_shared_counts, n_top_genes=n_top_genes)
    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
    return adata
def compute_velocity(adata, mode: str = 'stochastic'):
    """
    Compute RNA velocity using the specified mode.
    
    Args:
        adata: AnnData object prepared for velocity analysis.
        mode: Velocity mode ('deterministic', 'stochastic', or 'dynamical').
        
    Returns:
        AnnData object with computed velocity.
    """
    if mode == 'dynamical':
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode='dynamical')
    elif mode == 'deterministic':
        scv.tl.velocity(adata, mode='deterministic')
    else:
        scv.tl.velocity(adata, mode='stochastic')
    
    # Calculate velocity graph and pseudotime
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_pseudotime(adata)
    
    return adata
def plot_velocity_stream(adata, basis: str = 'umap', color: str = 'velocity_pseudotime', 
                        figdir: str = "figures", filename: str = "velocity_stream.png"):
    """
    Plot velocity stream plot.

    Args:
        adata: AnnData object with velocity results.
        basis: Basis for embedding.
        color: Variable to color by.
        figdir: Directory to save figures.
        filename: Name of the output file.
    """
    scv.pl.velocity_embedding_stream(adata, basis=basis, color=color, show=False, save=False)
    save_fig(None, filename, figdir, plt_figure=True)

def plot_velocity_arrows(adata, basis: str = 'umap', color: str = 'velocity_pseudotime', 
                        figdir: str = "figures", filename: str = "velocity_arrows.png"):
    """
    Plot velocity arrows on embedding.

    Args:
        adata: AnnData object with velocity results.
        basis: Basis for embedding.
        color: Variable to color by.
        figdir: Directory to save figures.
        filename: Name of the output file.
    """
    scv.pl.velocity_embedding(adata, basis=basis, color=color, show=False, save=False)
    save_fig(None, filename, figdir, plt_figure=True)
def plot_velocity_pseudotime(adata, basis: str = 'umap', 
                           figdir: str = "figures", filename: str = "velocity_pseudotime.png"):
    """
    Plot velocity pseudotime on embedding.

    Args:
        adata: AnnData object with velocity results.
        basis: Basis for embedding.
        figdir: Directory to save figures.
        filename: Name of the output file.
    """
    scv.pl.scatter(adata, basis=basis, color='velocity_pseudotime', 
                  show=False, save=False)
    save_fig(None, filename, figdir, plt_figure=True)

def identify_driver_genes(adata, n_genes: int = 30, 
                         figdir: str = "figures", filename: str = "driver_genes.png"):
    """
    Identify and plot top driver genes.

    Args:
        adata: AnnData object with velocity results.
        n_genes: Number of top genes to identify.
        figdir: Directory to save figures.
        filename: Name of the output file.
    """
    scv.tl.rank_velocity_genes(adata, n_genes=n_genes)
    df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
    
    scv.pl.velocity_genes(adata, show=False, save=False)
    save_fig(None, filename, figdir, plt_figure=True)
    
    return df
def plot_phase_portraits(adata, genes: List[str], 
                       figdir: str = "figures", ncols: int = 2):
    """
    Plot phase portraits for specified genes.

    Args:
        adata: AnnData object with velocity results.
        genes: List of genes to plot.
        figdir: Directory to save figures.
        ncols: Number of columns in the plot grid.
    """
    os.makedirs(figdir, exist_ok=True)
    
    # Filter genes to ensure they exist in the dataset
    valid_genes = [gene for gene in genes if gene in adata.var_names]
    
    if not valid_genes:
        print("None of the specified genes found in the dataset")
        return
    
    scv.pl.velocity(adata, valid_genes, ncols=ncols, show=False, save=False)
    save_fig(None, "phase_portraits.png", figdir, plt_figure=True)
