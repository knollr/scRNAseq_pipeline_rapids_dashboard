"""
components/clustering_umap.py
----------------------------------
Interactive UMAP and cluster statistics plot.

Works for both AnnData and MuData (via adata = mdata.mod["rna"]).
"""

import plotly.express as px
from dash import html, dcc, Input, Output
import pandas as pd
import numpy as np


def layout(adata):
    """Return layout for clustering + UMAP visualization."""

    # detect available clustering columns
    cluster_cols = [col for col in adata.obs.columns if col.startswith("leiden") or col.startswith("louvain")]
    if not cluster_cols:
        cluster_cols = ["no_cluster_found"]

    return html.Div(
        [
            html.H3("Cluster UMAP Visualization"),
            dcc.Dropdown(
                id="cluster-dropdown",
                options=[{"label": c, "value": c} for c in cluster_cols],
                value=cluster_cols[0],
                clearable=False,
                style={"width": "50%"},
            ),
            html.Div(
                [
                    dcc.Graph(id="umap-plot"),
                    dcc.Graph(id="cluster-barplot"),
                ],
                style={"display": "grid", "gridTemplateColumns": "1fr 1fr", "gap": "1rem"},
            ),
        ]
    )


def register_callbacks(app, adata):
    """Register all callbacks for interactive behavior."""

    @app.callback(
        [Output("umap-plot", "figure"), Output("cluster-barplot", "figure")],
        Input("cluster-dropdown", "value"),
    )
    def update_plots(cluster_col):
        if cluster_col not in adata.obs.columns:
            return px.scatter(title="No clustering column found."), px.bar(title="No data")

        df = adata.obs.copy()
        df["UMAP1"] = adata.obsm["X_umap"][:, 0]
        df["UMAP2"] = adata.obsm["X_umap"][:, 1]

        df_sample = df.sample(n=min(10000, df.shape[0]), random_state=0)

        # --- UMAP scatter ---
        fig_umap = px.scatter(
            df_sample,
            x="UMAP1",
            y="UMAP2",
            color=cluster_col,
            title=f"UMAP colored by {cluster_col}",
            opacity=0.7,
            width=600,
            height=600,
            render_mode="webgl",  # 🚀 GPU acceleration
        )
        fig_umap.update_traces(marker=dict(size=3, line=dict(width=0)))

        # --- Cluster counts ---
        cluster_counts = df[cluster_col].value_counts().sort_index().reset_index()
        cluster_counts.columns = ["Cluster", "Cell count"]

        fig_bar = px.bar(
            cluster_counts,
            x="Cluster",
            y="Cell count",
            title=f"Cells per {cluster_col}",
            text_auto=True,
            width=600,
            height=600,
        )

        return fig_umap, fig_bar
