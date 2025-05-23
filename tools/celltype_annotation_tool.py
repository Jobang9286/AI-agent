# tools/celltype_annotation_tool.py

from server import mcp
import scanpy as sc
import pandas as pd
import numpy as np
from utils.plot_utils import save_fig

@mcp.tool()
def annotate_cell_types(
    adata_path: str, 
    reference_markers: str = None,
    method: str = "manual", 
    cluster_key: str = "leiden",
    organism: str = "human",
    figdir: str = "figures"
) -> str:
    """
    Annotate cell types based on marker genes.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        reference_markers: Path to CSV file with reference markers (gene,cell_type).
        method: Annotation method ('manual', 'automated').
        cluster_key: Key in adata.obs containing cluster information.
        organism: Organism ('human' or 'mouse').
        figdir: Directory to save annotation visualizations.

    Returns:
        Path to the .h5ad file with cell type annotations.
    """
    adata = sc.read_h5ad(adata_path)
    
    if method == "manual" and reference_markers:
        # Load marker genes from reference file
        markers_df = pd.read_csv(reference_markers)
        
        # Create a dictionary mapping cell types to marker genes
        cell_type_markers = {}
        for _, row in markers_df.iterrows():
            if row['cell_type'] not in cell_type_markers:
                cell_type_markers[row['cell_type']] = []
            cell_type_markers[row['cell_type']].append(row['gene'])
        
        # Score cells based on the expression of marker genes
        for cell_type, markers in cell_type_markers.items():
            # Filter for markers that exist in the dataset
            markers = [m for m in markers if m in adata.var_names]
            if markers:
                sc.tl.score_genes(adata, markers, score_name=f"{cell_type}_score")
        
        # Create a plot for each cell type score
        score_columns = [c for c in adata.obs.columns if c.endswith('_score')]
        if score_columns and 'X_umap' in adata.obsm:
            sc.pl.umap(adata, color=score_columns, show=False, save="_cell_type_scores.png")
        
        # Add a categorical cell type annotation based on the highest score
        if score_columns:
            adata.obs['cell_type'] = 'Unknown'
            for idx, row in adata.obs.iterrows():
                scores = [(score_col.replace('_score', ''), row[score_col]) for score_col in score_columns]
                best_type, _ = max(scores, key=lambda x: x[1])
                adata.obs.at[idx, 'cell_type'] = best_type
            
            # Plot the annotated cell types
            if 'X_umap' in adata.obsm:
                sc.pl.umap(adata, color='cell_type', show=False, save="_cell_type_annotation.png")
                
    elif method == "automated" or not reference_markers:
        # Simple automated annotation based on common brain cell markers
        if organism == "mouse":
            brain_markers = {
                'Neurons': ['Snap25', 'Syn1', 'Rbfox3', 'Tubb3'],
                'Astrocytes': ['Gfap', 'Slc1a3', 'Aqp4', 'Aldoc'],
                'Oligodendrocytes': ['Mbp', 'Olig2', 'Sox10', 'Cnp'],
                'Microglia': ['Cx3cr1', 'P2ry12', 'Csf1r', 'Hexb'],
                'Endothelial': ['Cldn5', 'Pecam1', 'Flt1'],
                'Neural_progenitors': ['Sox2', 'Nestin', 'Pax6', 'Sox4']
            }
        else:  # human
            brain_markers = {
                'Neurons': ['SNAP25', 'SYN1', 'RBFOX3', 'TUBB3'],
                'Astrocytes': ['GFAP', 'SLC1A3', 'AQP4', 'ALDOC'],
                'Oligodendrocytes': ['MBP', 'OLIG2', 'SOX10', 'CNP'],
                'Microglia': ['CX3CR1', 'P2RY12', 'CSF1R', 'HEXB'],
                'Endothelial': ['CLDN5', 'PECAM1', 'FLT1'],
                'Neural_progenitors': ['SOX2', 'NES', 'PAX6', 'SOX4']
            }
        
        # Score cells for each cell type
        for cell_type, markers in brain_markers.items():
            markers = [m for m in markers if m in adata.var_names]
            if markers:
                sc.tl.score_genes(adata, markers, score_name=f"{cell_type}_score")
        
        # Assign cell types based on highest score
        score_columns = [c for c in adata.obs.columns if c.endswith('_score')]
        if score_columns:
            adata.obs['cell_type'] = 'Unknown'
            for idx, row in adata.obs.iterrows():
                scores = [(score_col.replace('_score', ''), row[score_col]) for score_col in score_columns if not np.isnan(row[score_col])]
                if scores:
                    best_type, _ = max(scores, key=lambda x: x[1])
                    adata.obs.at[idx, 'cell_type'] = best_type
            
            # Plot results if UMAP exists
            if 'X_umap' in adata.obsm:
                sc.pl.umap(adata, color='cell_type', show=False, save="_automated_annotation.png")
    
    # Save the annotated object
    out_path = adata_path.replace(".h5ad", "_annotated.h5ad")
    adata.write(out_path)
    return out_path