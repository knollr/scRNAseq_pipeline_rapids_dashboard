#!/usr/bin/env python3
"""
app.py — Unified dashboard for single-cell pipeline visualization.

Run inside Singularity container:
    singularity exec -B $(pwd):/pipeline --nv scrnaseq_pipeline_latest.sif \
        python /pipeline/dashboard/app.py

Assumes Snakemake created the symlink:
    results/latest/merged_annotated.(h5ad|h5mu)
"""
import components.silence_warnings_app
import os
import sys
import dash
from dash import html, dcc
import scanpy as sc
import muon as mu
import pandas as pd
import argparse

# ---------------------------------------------------------------------
# Craete DataWrapper class as simple container for AnnData and MuData
# ---------------------------------------------------------------------
class DataWrapper:
    def __init__(self, adata, mdata=None):
        self.adata = adata       # the AnnData to use for all RNA plots
        self.mdata = mdata       # full MuData object if available
        self.is_mdata = mdata is not None



# ---------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Single-cell dashboard")
parser.add_argument(
    "--results-dir",
    type=str,
    default=None,
    help="Path to results directory (defaults to results/latest)",
)
args = parser.parse_args()

# ---------------------------------------------------------------------
# Resolve results directory
# ---------------------------------------------------------------------
DEFAULT_RESULTS_DIR = "/pipeline/results/latest"

RESULTS_DIR = (
    args.results_dir
    if args.results_dir is not None
    else DEFAULT_RESULTS_DIR
)

print(f"📂 Using results directory: {RESULTS_DIR}")

# ---------------------------------------------------------------------
# Load annotated dataset
# ---------------------------------------------------------------------
ANNOT_PATH_H5AD = os.path.join(RESULTS_DIR, "merged_annotated.h5ad")
ANNOT_PATH_H5MU = os.path.join(RESULTS_DIR, "merged_annotated.h5mu")
ANNOT_PATH_CSV = os.path.join(RESULTS_DIR, "doublet_qc_dataframe_pre_filtering.csv")

adata = None
if os.path.exists(ANNOT_PATH_H5MU):
    print(f"📂 Loading MuData from {ANNOT_PATH_H5MU}")
    mdata = mu.read(ANNOT_PATH_H5MU)
    adata = mdata.mod["rna"]  
    data = DataWrapper(adata=adata, mdata=mdata) # wrap both
elif os.path.exists(ANNOT_PATH_H5AD):
    print(f"📂 Loading AnnData from {ANNOT_PATH_H5AD}")
    adata = sc.read_h5ad(ANNOT_PATH_H5AD)
    data = DataWrapper(adata=adata, mdata=None)  # only adata
else:
    sys.exit(f"❌ No merged_annotated.h5ad or .h5mu found in {RESULTS_DIR}")
    
doublet_df = pd.read_csv(
    ANNOT_PATH_CSV,
    index_col="cell_id"
)

data.doublet_df = doublet_df   # attach to your DataLoader object    

print(f"✅ Loaded dataset with {data.adata.n_obs} cells × {data.adata.n_vars} genes.")


# ---------------------------------------------------------------------
# Import dashboard components
# ---------------------------------------------------------------------
from components import (
    qc_plots,
    clustering_umap,
    celltype_summary,
)

# ---------------------------------------------------------------------
# Initialize Dash app
# ---------------------------------------------------------------------
app = dash.Dash(
    __name__,
    title="Single-cell Dashboard",
    suppress_callback_exceptions=True,
)
server = app.server

# ---------------------------------------------------------------------
# App layout
# ---------------------------------------------------------------------
app.layout = html.Div(
    [
        html.H1("🧬 Single-cell Analysis Dashboard", style={"textAlign": "center"}),
        html.Hr(),
        dcc.Tabs(
            id="tabs",
            value="tab-qc",
            children=[
                dcc.Tab(label="🔍 QC Overview", value="tab-qc"),
                dcc.Tab(label="🌐 Clustering & UMAP", value="tab-cluster"),
                dcc.Tab(label="🧫 Cell Type Summary", value="tab-celltype"),
            ],
        ),
        html.Div(id="tab-content", style={"margin": "2rem"}),
    ]
)

# ---------------------------------------------------------------------
# Callbacks for tab content
# ---------------------------------------------------------------------
@app.callback(
    dash.Output("tab-content", "children"),
    dash.Input("tabs", "value"),
)
def render_tab(tab):
    if tab == "tab-qc":
        return qc_plots.layout(data)
    elif tab == "tab-cluster":
        return clustering_umap.layout(data)
    elif tab == "tab-celltype":
        return celltype_summary.layout(data)
    else:
        return html.Div("Select a tab to view content.")

# ---------------------------------------------------------------------
# Register component callbacks
# ---------------------------------------------------------------------
qc_plots.register_callbacks(app, data)
clustering_umap.register_callbacks(app, data)
celltype_summary.register_callbacks(app, data)

# ---------------------------------------------------------------------
# Run dashboard
# ---------------------------------------------------------------------
if __name__ == "__main__":
    print("🚀 Dashboard running at http://0.0.0.0:8050")
    app.run(host="0.0.0.0", port=8050, debug=False)
