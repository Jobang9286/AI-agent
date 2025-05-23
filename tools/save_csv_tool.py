# tools/save_csv_tool.py

from server import mcp
import scanpy as sc
import os

@mcp.tool()
def save_csv_tables(
    adata_path: str,
    outdir: str = "exports",
    save_obs: bool = True,
    save_var: bool = False,
    save_ranked_genes: bool = True
) -> str:
    """
    Export key AnnData tables (obs, var, ranked_genes) to CSV files.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        outdir: Directory to save CSV files.
        save_obs: Whether to save adata.obs.
        save_var: Whether to save adata.var.
        save_ranked_genes: Whether to save rank_genes_groups as a summary table.

    Returns:
        Path to the output directory containing exported CSVs.
    """
    adata = sc.read_h5ad(adata_path)
    os.makedirs(outdir, exist_ok=True)

    if save_obs:
        adata.obs.to_csv(os.path.join(outdir, "obs.csv"))

    if save_var:
        adata.var.to_csv(os.path.join(outdir, "var.csv"))

    if save_ranked_genes and "rank_genes_groups" in adata.uns:
        ranked = sc.get.rank_genes_groups_df(adata, group=None)
        ranked.to_csv(os.path.join(outdir, "ranked_genes.csv"), index=False)

    return outdir
