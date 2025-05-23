from server import mcp

import scanpy as sc

from utils.prep_utils import qc_violin
from utils.prep_utils import filtering


@mcp.tool()
def qc_plot(file_path: str):
  """
  A violin plot of some of the computed quality measures:
    1. the number of genes expressed in the count matrix
    2. the total counts per cell
    3. the percentage of counts in mitochondrial genes

  Args:
    file_path

  Returns:
    None
  """
  adata = sc.read_h5ad(file_path)
  qc_violin(adata)
    
  return None


@mcp.tool()
def prep_filtering(file_path: str, n_genes_by_counts: int, pct_counts_mt: int) -> str:
  """
  This is a basic sc data filtering based on the 'n_genes_by_counts' and 'pct_counts_mt'
  Remove cells that have too many the number of genes expressed in the count matrix (i.e., n_genes_by_counts) or
  too many mitochondrial genes expressed (i.e., pct_counts_mt)

  Args:
    AnnData object
    n_genes_by_counts: the number of genes expressed in the count matrix
    pct_counts_mt: the percentage of counts in mitochondrial genes

  Returns:
    AnnData object
  """
  adata = sc.read_h5ad(file_path)
  adata = filtering(adata, n_genes_by_counts, pct_counts_mt)

  # 기본 정보 출력
  print("데이터 모양:", adata.shape)
  print(f"관측값 수: {adata.n_obs}, 변수 수: {adata.n_vars}")

  savedir2 = file_path.replace('.h5ad', '_filtered.h5ad')
  adata.write_h5ad(savedir2)

  return savedir2