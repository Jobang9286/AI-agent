
from server import mcp  
import scanpy as sc

@mcp.tool()
def normalize_data(adata_path: str, target_sum: float = 1e4) -> str:
    """
    Normalize and log1p-transform the data.

    Args:
        adata_path: Path to the AnnData .h5ad file.
        target_sum: Total counts per cell after normalization.

    Returns:
        Path to the updated AnnData object.
    """
    adata = sc.read_h5ad(adata_path)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    out_path = adata_path.replace(".h5ad", "_normalized.h5ad")
    adata.write(out_path)
    return out_path
