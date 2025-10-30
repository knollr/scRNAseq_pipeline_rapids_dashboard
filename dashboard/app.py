#!/usr/bin/env python3
"""
app.py — Entry point for the interactive dashboard.

Automatically detects and loads the latest single-cell output (AnnData or MuData)
from the Snakemake pipeline results.

Usage (inside Singularity container):
    singularity exec -B $(pwd):/pipeline --nv scrnaseq_pipeline_latest.sif \
        python /pipeline/dashboard/app.py
"""

import os
import sys
import dash
from dash import html, dcc
import scanpy as sc
import muon as mu

# ---------------------------------------------------------------------
# Detect and load latest dataset
# ---------------------------------------------------------------------
RESULTS_DIR = "/pipeline/results/latest"
if not os.path.exists(RESULTS_DIR):
    sys.exit(f"❌ Error: Expected symlink {RESULTS_DIR} not found. Did the pipeline finish?")

# Detect file type (priority: .h5mu > .h5ad)
mdata_path = os.path.join(RESULTS_DIR, "merged_annotated.h5mu")
adata_path = os.path.join(RESULTS_DIR, "merged_annotated.h5ad")

adata = None
mdata = None

if os.path.exists(mdata_path):
    print(f"📂 Loading MuData: {mdata_path}")
    mdata = mu.read(mdata_path)
    # assume RNA modality exists
    adata = mdata.mod["rna"]
elif os.path.exists(adata_path):
    print(f"📂 Loading AnnData: {adata_path}")
    adata = sc.read_h5ad(adata_path)
else:
    sys.exit(f"❌ No annotated file found in {RESULTS_DIR} (expected .h5ad or .h5mu)")

print(f"✅ Loaded object with {adata.n_obs} cells × {adata.n_vars} genes.")

# ---------------------------------------------------------------------
# Import components
# ---------------------------------------------------------------------
from components import clustering_umap  # imports layout + callbacks

# ---------------------------------------------------------------------
# Initialize Dash app
# ---------------------------------------------------------------------
app = dash.Dash(__name__, title="Single-cell Dashboard", suppress_callback_exceptions=True)
server = app.server  # for deployment use

app.layout = html.Div(
    [
        html.H2("🧬 Single-cell Analysis Dashboard", style={"textAlign": "center"}),
        html.Hr(),
        clustering_umap.layout(adata),
    ],
    style={"margin": "2rem"},
)

# Register callbacks after layout
clustering_umap.register_callbacks(app, adata)

# ---------------------------------------------------------------------
# Run server
# ---------------------------------------------------------------------
if __name__ == "__main__":
    app.run(debug=False, host="0.0.0.0", port=8050)
