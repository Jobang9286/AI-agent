# utils/pathway_utils.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from typing import List, Dict, Optional, Union, Tuple
import scanpy as sc
from utils.plot_utils import save_fig

def run_enrichr(gene_list: List[str], gene_set: str = "GO_Biological_Process_2021", 
               organism: str = "Human", description: str = "Enrichment"):
    """
    Run Enrichr enrichment analysis on a gene list.

    Args:
        gene_list: List of gene symbols.
        gene_set: Gene set library to use.
        organism: 'Human' or 'Mouse'.
        description: Description for the analysis.

    Returns:
        Enrichment results dataframe.
    """
    try:
        import gseapy as gp
        enr = gp.enrichr(gene_list=gene_list,
                         gene_sets=gene_set,
                         organism=organism,
                         description=description,
                         outdir=None,
                         no_plot=True)
        
        return enr.results
    except ImportError:
        print("gseapy package not installed. Please install it with 'pip install gseapy'")
        return pd.DataFrame()
    except Exception as e:
        print(f"Error running enrichment analysis: {e}")
        return pd.DataFrame()
def plot_enrichment_barplot(enrichment_results: pd.DataFrame, n_top: int = 10, 
                           figdir: str = "figures", filename: str = "enrichment_barplot.png",
                           title: str = "Top Enriched Pathways"):
    """
    Create a bar plot of top enriched pathways.

    Args:
        enrichment_results: Enrichment results dataframe from run_enrichr.
        n_top: Number of top pathways to include.
        figdir: Directory to save figure.
        filename: Filename for the saved figure.
        title: Title for the plot.
    """
    # Filter and sort results
    filtered_results = enrichment_results[enrichment_results['Adjusted P-value'] < 0.05]
    filtered_results = filtered_results.sort_values('Adjusted P-value').head(n_top)
    
    if len(filtered_results) == 0:
        print("No significant results found.")
        return
    
    # Create bar plot
    plt.figure(figsize=(10, 6))
    plt.barh(range(len(filtered_results)), 
            -np.log10(filtered_results['Adjusted P-value']), 
            color='skyblue')
    plt.yticks(range(len(filtered_results)), filtered_results['Term'])
    plt.xlabel('-log10(Adjusted P-value)')
    plt.title(title)
    plt.tight_layout()
    
    os.makedirs(figdir, exist_ok=True)
    plt.savefig(f"{figdir}/{filename}", dpi=300, bbox_inches='tight')
    plt.close()
def get_marker_genes(adata, group: str, groupby: str = 'leiden', 
                    adj_pval_threshold: float = 0.05, 
                    log2fc_threshold: float = 0.25) -> List[str]:
    """
    Extract significant marker genes for a specific cluster.
    
    Args:
        adata: AnnData object with rank_genes_groups results
        group: The group/cluster to get markers for
        groupby: The key used for grouping in rank_genes_groups
        adj_pval_threshold: Adjusted p-value threshold for significance
        log2fc_threshold: Log2 fold change threshold
        
    Returns:
        List of significant marker genes
    """
    if 'rank_genes_groups' not in adata.uns:
        print("No rank_genes_groups found in adata. Running wilcoxon test...")
        sc.tl.rank_genes_groups(adata, groupby=groupby, method='wilcoxon')
    
    # Get results dataframe for the specified group
    markers_df = sc.get.rank_genes_groups_df(adata, group=group)
    
    # Filter by significance and fold change
    significant_markers = markers_df[
        (markers_df['pvals_adj'] < adj_pval_threshold) & 
        (markers_df['logfoldchanges'] > log2fc_threshold)
    ]
    
    return significant_markers['names'].tolist()
def pathway_summary(enrichment_results_dict: Dict[str, pd.DataFrame], 
                   top_n: int = 3, 
                   figdir: str = "figures", 
                   filename: str = "pathway_summary.csv"):
    """
    Create a summary of top pathways for each group.
    
    Args:
        enrichment_results_dict: Dictionary mapping groups to enrichment results
        top_n: Number of top pathways to include per group
        figdir: Directory to save results
        filename: Name of the output file
        
    Returns:
        DataFrame with summary information
    """
    summary_rows = []
    
    for group, results in enrichment_results_dict.items():
        if len(results) == 0 or 'Adjusted P-value' not in results.columns:
            continue
            
        # Filter and sort results
        filtered_results = results[results['Adjusted P-value'] < 0.05]
        filtered_results = filtered_results.sort_values('Adjusted P-value')
        
        # Get top pathways
        for i, (_, pathway) in enumerate(filtered_results.head(top_n).iterrows()):
            summary_rows.append({
                'Group': group,
                'Rank': i + 1,
                'Pathway': pathway['Term'],
                'Adjusted P-value': pathway['Adjusted P-value'],
                'Genes': pathway['Genes']
            })
    
    # Create dataframe
    summary_df = pd.DataFrame(summary_rows)
    
    if len(summary_df) > 0:
        # Save to CSV
        os.makedirs(figdir, exist_ok=True)
        summary_df.to_csv(f"{figdir}/{filename}", index=False)
    
    return summary_df
def create_heatmap_from_pathways(adata, pathways_dict: Dict[str, List[str]], 
                               groupby: str = 'leiden', 
                               figdir: str = "figures", 
                               filename: str = "pathway_heatmap.png"):
    """
    Create a heatmap showing pathway activity across groups.
    
    Args:
        adata: AnnData object
        pathways_dict: Dictionary mapping pathway names to lists of genes
        groupby: Column in adata.obs to group cells by
        figdir: Directory to save figure
        filename: Name of output file
    """
    # Create a dataframe to store pathway scores
    pathway_scores = pd.DataFrame(index=adata.obs[groupby].unique())
    
    # Calculate score for each pathway
    for pathway_name, gene_list in pathways_dict.items():
        # Filter to genes present in the dataset
        valid_genes = [gene for gene in gene_list if gene in adata.var_names]
        
        if len(valid_genes) == 0:
            print(f"No genes from pathway '{pathway_name}' found in dataset")
            continue
        
        # Score cells by pathway activity
        temp_adata = adata.copy()
        sc.tl.score_genes(temp_adata, valid_genes, score_name='pathway_score')
        
        # Calculate mean score per group
        group_means = temp_adata.obs.groupby(groupby)['pathway_score'].mean()
        pathway_scores[pathway_name] = group_means
    
    if pathway_scores.empty or pathway_scores.shape[1] == 0:
        print("No valid pathways to plot")
        return
    
    # Create heatmap
    plt.figure(figsize=(10, 8))
    im = plt.imshow(pathway_scores, cmap='viridis')
    plt.colorbar(im, label='Mean Pathway Score')
    
    # Add labels
    plt.yticks(range(len(pathway_scores.index)), pathway_scores.index)
    plt.xticks(range(len(pathway_scores.columns)), pathway_scores.columns, rotation=45, ha='right')
    
    plt.title('Pathway Activity Across Cell Groups')
    plt.tight_layout()
    
    # Save figure
    os.makedirs(figdir, exist_ok=True)
    plt.savefig(f"{figdir}/{filename}", dpi=300, bbox_inches='tight')
    plt.close()
