# tools/scale_tool.py

from server import mcp
import scanpy as sc

@mcp.tool()
def scale_data(adata_path: str, max_value: float = 10.0) -> str:
    """
    Scale each gene to unit variance and zero mean. Clip values above max_value.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        max_value: Value to clip extreme values after scaling.

    Returns:
        Path to the scaled .h5ad file.
    """
    adata = sc.read_h5ad(adata_path)
    sc.pp.scale(adata, max_value=max_value)

    out_path = adata_path.replace(".h5ad", "_scaled.h5ad")
    adata.write(out_path)
    return out_path
