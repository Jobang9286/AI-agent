# tools/hvg_tool.py

from server import mcp
import scanpy as sc

@mcp.tool()
def select_highly_variable_genes(adata_path: str, min_mean: float = 0.0125, max_mean: float = 3, min_disp: float = 0.5) -> str:
    """
    Select highly variable genes for downstream analysis.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        min_mean: Minimum mean expression threshold.
        max_mean: Maximum mean expression threshold.
        min_disp: Minimum dispersion threshold.

    Returns:
        Path to the .h5ad file with highly variable genes annotated and filtered.
    """
    adata = sc.read_h5ad(adata_path)
    sc.pp.highly_variable_genes(adata, min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)
    adata = adata[:, adata.var.highly_variable]

    out_path = adata_path.replace(".h5ad", "_hvg.h5ad")
    adata.write(out_path)
    return out_path
