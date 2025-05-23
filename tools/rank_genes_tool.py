# tools/rank_genes_tool.py

from server import mcp
import scanpy as sc
from utils.plot_utils import save_fig

@mcp.tool()
def find_marker_genes(
    adata_path: str,
    groupby: str = "leiden",
    method: str = "t-test",
    figdir: str = "figures"
) -> str:
    """
    Identify marker genes for each cluster and save a rank_genes_groups plot.

    Args:
        adata_path: Path to the input AnnData (.h5ad) file.
        groupby: Column in adata.obs to group cells for comparison.
        method: Statistical test method (e.g., 't-test', 'wilcoxon').
        figdir: Directory to save the marker plot.

    Returns:
        Path to the .h5ad file with marker gene results.
    """
    import subprocess
    import sys
    import os
    
    # Create a standalone script that runs in the correct environment
    script_content = f'''
import scanpy as sc
import matplotlib.pyplot as plt
import warnings
import os
warnings.filterwarnings("ignore")

# Load data and compute marker genes
adata = sc.read("{adata_path}")
sc.tl.rank_genes_groups(adata, groupby="{groupby}", method="{method}")

# Save plot if figdir provided
if "{figdir}":
    os.makedirs("{figdir}", exist_ok=True)
    try:
        sc.pl.rank_genes_groups(adata, sharey=False, show=False)
        fig = plt.gcf()
        fig.savefig(os.path.join("{figdir}", "rank_genes_{groupby}.png"), dpi=300, bbox_inches="tight")
        plt.close()
        print("Plot saved successfully")
    except Exception as e:
        print(f"Plot save failed: {{e}}")

# Save data
out_path = "{adata_path}".replace(".h5ad", "_ranked.h5ad")
adata.write(out_path)
print(f"Data saved to {{out_path}}")
'''
    
    # Write and execute the script
    script_path = "/tmp/rank_genes_script.py"
    with open(script_path, "w") as f:
        f.write(script_content)
    
    try:
        # Execute using the virtual environment python
        venv_python = "/Users/bang/Desktop/scanpy_bg/.venv/bin/python"
        if os.path.exists(venv_python):
            result = subprocess.run([venv_python, script_path], 
                                  capture_output=True, text=True, timeout=60)
        else:
            result = subprocess.run([sys.executable, script_path], 
                                  capture_output=True, text=True, timeout=60)
        
        print("Script output:", result.stdout)
        if result.stderr:
            print("Script errors:", result.stderr)
            
        if result.returncode != 0:
            raise Exception(f"Script failed with return code {result.returncode}")
    
    finally:
        # Clean up
        if os.path.exists(script_path):
            os.remove(script_path)
    
    out_path = adata_path.replace(".h5ad", "_ranked.h5ad")
    return out_path
