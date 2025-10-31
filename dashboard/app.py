#!/usr/bin/env python3
"""
app.py — Unified dashboard for single-cell pipeline visualization.

Run inside Singularity container:
    singularity exec -B $(pwd):/pipeline --nv scrnaseq_pipeline_latest.sif \
        python /pipeline/dashboard/app.py

Assumes Snakemake created the symlink:
    results/latest/merged_annotated.(h5ad|h5mu)
"""

import os
import sys
import dash
from dash import html, dcc
import scanpy as sc
import muon as mu

# ---------------------------------------------------------------------
# Load annotated dataset
# ---------------------------------------------------------------------
RESULTS_DIR = "/pipeline/results/latest"
ANNOT_PATH_H5AD = os.path.join(RESULTS_DIR, "merged_annotated.h5ad")
ANNOT_PATH_H5MU = os.path.join(RESULTS_DIR, "merged_annotated.h5mu")

adata = None
if os.path.exists(ANNOT_PATH_H5MU):
    print(f"📂 Loading MuData from {ANNOT_PATH_H5MU}")
    mdata = mu.read(ANNOT_PATH_H5MU)
    adata = mdata.mod["rna"]  # Focus on RNA modality
elif os.path.exists(ANNOT_PATH_H5AD):
    print(f"📂 Loading AnnData from {ANNOT_PATH_H5AD}")
    adata = sc.read_h5ad(ANNOT_PATH_H5AD)
else:
    sys.exit(f"❌ No merged_annotated.h5ad or .h5mu found in {RESULTS_DIR}")

print(f"✅ Loaded dataset with {adata.n_obs} cells × {adata.n_vars} genes.")

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
        return qc_plots.layout(adata)
    elif tab == "tab-cluster":
        return clustering_umap.layout(adata)
    elif tab == "tab-celltype":
        return celltype_summary.layout(adata)
    else:
        return html.Div("Select a tab to view content.")

# ---------------------------------------------------------------------
# Register component callbacks
# ---------------------------------------------------------------------
qc_plots.register_callbacks(app, adata)
clustering_umap.register_callbacks(app, adata)
celltype_summary.register_callbacks(app, adata)

# ---------------------------------------------------------------------
# Run dashboard
# ---------------------------------------------------------------------
if __name__ == "__main__":
    print("🚀 Dashboard running at http://0.0.0.0:8050")
    app.run(host="0.0.0.0", port=8050, debug=False)
