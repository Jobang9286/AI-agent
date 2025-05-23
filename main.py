from server import mcp

from tools.prep_tool import qc_plot, prep_filtering
from tools.read_tool import read_sc
from tools.normalized_tool import normalize_data
from tools.hvg_tool import select_highly_variable_genes
from tools.scale_tool import scale_data
from tools.pca_tool import compute_pca
from tools.neighbors_tool import compute_neighbors
from tools.clustering_tool import leiden_clustering
from tools.umap_tool import compute_umap
from tools.rank_genes_tool import find_marker_genes
import tools.dotplot_tool
import tools.heatmap_tool
import tools.save_csv_tool

from tools.celltype_annotation_tool import annotate_cell_types
from tools.pseudotime_tool import compute_pseudotime
from tools.velocity_tool import compute_rna_velocity
from tools.gene_module_tool import analyze_gene_modules
from tools.scanpy_pathway_tool import pathway_analysis
from tools.bbknn_tool import compute_bbknn
from tools.harmony_tool import compute_harmony

if __name__ == "__main__":
    mcp.run()
