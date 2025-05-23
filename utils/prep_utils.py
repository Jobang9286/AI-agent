import scanpy as sc
import anndata as ad


def read_10x_mtx(filename: str) -> ad.AnnData:

  adata = sc.read_10x_mtx(
      filename,  # the directory with the `.mtx` file
      var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
  )
  print("데이터 모양:", adata.shape)
  print(f"관측값 수: {adata.n_obs}, 변수 수: {adata.n_vars}")

  adata.var_names_make_unique()
  # For mouse data, use "mt-" instead of "MT-" (human)
  adata.var["mt"] = adata.var_names.str.startswith("mt-")
  sc.pp.calculate_qc_metrics(
      adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
  )
  
  return adata


def qc_violin(adata: ad.AnnData):

  print("plot vioin...")
  sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save="_qc.png"
  )

  return None


def filtering(adata: ad.AnnData, n_genes_by_counts: int=2500, pct_counts_mt: int=5) -> ad.AnnData:

  adata = adata[adata.obs.n_genes_by_counts < n_genes_by_counts, :]
  adata = adata[adata.obs.pct_counts_mt < pct_counts_mt, :].copy()

  return adata