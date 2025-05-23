from server import mcp

import scanpy as sc

from utils.prep_utils import read_10x_mtx


@mcp.tool()
def read_sc(
  filename: str, 
  savedir: str
) -> str:
  """
  Read in the count matrix (i.e., matrix.mtx) into an AnnData object (https://anndata.readthedocs.io/en/stable/), 
  which holds many slots for annotations and different representations of the data.
  It also comes with its own HDF5-based file format: .h5ad.

  Args:
    filename: the directory with the `.mtx` file

  Returns:
    AnnData object
  """
  adata = read_10x_mtx(filename)
  
  savedir = f"{savedir}/adata.h5ad"
  adata.write_h5ad(savedir)
  
  return savedir