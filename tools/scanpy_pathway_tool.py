# tools/scanpy_pathway_tool.py

from server import mcp
import scanpy as sc
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from utils.plot_utils import save_fig

@mcp.tool()
def pathway_analysis(
    adata_path: str,
    groupby: str = "leiden",
    gene_set: str = "GO_Biological_Process_2021",
    organism: str = "Human",
    figdir: str = "figures"
) -> str:
    """
    Perform pathway and gene set enrichment analysis on marker genes.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        groupby: Column in adata.obs to group cells for comparison.
        gene_set: Gene set to use for enrichment (e.g., 'GO_Biological_Process_2021', 'KEGG_2021_Human').
        organism: Organism for gene sets ('Human' or 'Mouse').
        figdir: Directory to save pathway analysis results.

    Returns:
        Path to the .h5ad file with pathway analysis results.
    """
    adata = sc.read_h5ad(adata_path)
    
    # Check if rank_genes_groups has been run
    if 'rank_genes_groups' not in adata.uns:
        sc.tl.rank_genes_groups(adata, groupby=groupby, method='wilcoxon')
    
    # Create directory for saving results
    os.makedirs(figdir, exist_ok=True)
    
    # Store the enrichment results in a dict
    enrichment_results = {}
    
    # Get unique groups
    groups = adata.obs[groupby].unique()
    
    # Run enrichment analysis for each group
    for group in groups[:min(len(groups), 5)]:  # Limit to first 5 groups to avoid long runtime
        try:
            # Get marker genes for this group
            markers_df = sc.get.rank_genes_groups_df(adata, group=group)
            markers_df = markers_df[markers_df['pvals_adj'] < 0.05]  # Filter by adjusted p-value
            
            # If no significant genes, skip
            if len(markers_df) == 0:
                continue
            
            # Get gene list with log2 fold changes
            gene_list = markers_df.set_index('names')['logfoldchanges'].to_dict()
            
            # Run enrichment analysis
            import gseapy as gp
            enr = gp.enrichr(gene_list=list(gene_list.keys()),
                            gene_sets=gene_set,
                            organism=organism,
                            description=f'Cluster_{group}',
                            outdir=None,
                            no_plot=True)
            
            # Filter results by adjusted p-value
            enr_results = enr.results
            enr_results = enr_results[enr_results['Adjusted P-value'] < 0.05]
            
            # Save to dictionary
            enrichment_results[group] = enr_results            
            # Save top pathways to CSV
            enr_results.to_csv(f"{figdir}/enrichment_cluster_{group}.csv", index=False)
            
            # Plot top 10 pathways for this group
            if len(enr_results) > 0:
                # Sort by adjusted p-value
                enr_results = enr_results.sort_values('Adjusted P-value')
                
                # Take top 10 pathways
                top_pathways = enr_results.head(10)
                
                # Create bar plot
                plt.figure(figsize=(10, 6))
                plt.barh(range(len(top_pathways)), -np.log10(top_pathways['Adjusted P-value']), color='skyblue')
                plt.yticks(range(len(top_pathways)), top_pathways['Term'])
                plt.xlabel('-log10(Adjusted P-value)')
                plt.title(f'Top Pathways for Cluster {group}')
                plt.tight_layout()
                plt.savefig(f"{figdir}/enrichment_cluster_{group}_barplot.png", dpi=300, bbox_inches='tight')
                plt.close()
        
        except Exception as e:
            print(f"Error in enrichment analysis for group {group}: {e}")
    
    # Create a summary dataframe of top pathways per cluster
    summary_rows = []
    for group, results in enrichment_results.items():
        if len(results) > 0:
            top_path = results.iloc[0]
            summary_rows.append({
                'cluster': group,
                'top_pathway': top_path['Term'],
                'p_value': top_path['Adjusted P-value'],
                'genes': top_path['Genes']
            })    
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv(f"{figdir}/enrichment_summary.csv", index=False)
    
    # Store the enrichment results paths in adata.uns
    adata.uns['pathway_analysis'] = {
        'groupby': groupby,
        'gene_set': gene_set,
        'organism': organism,
        'result_files': [f for f in os.listdir(figdir) if f.startswith('enrichment_')]
    }
    
    # Save results
    out_path = adata_path.replace(".h5ad", "_pathway.h5ad")
    adata.write(out_path)
    return out_path