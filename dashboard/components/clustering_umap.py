# dashboard/components/clustering_umap.py

import plotly.express as px
from dash import html, dcc
import pandas as pd

def clustering_umap_component(adata):
    """
    Returns a Dash HTML Div containing an interactive UMAP scatter plot
    and a cluster count bar plot, with a dropdown for choosing cluster resolution.
    """

    # Extract available leiden columns
    cluster_cols = [col for col in adata.obs.columns if col.startswith("leiden_")]
    default_cluster = cluster_cols[0] if cluster_cols else None

    # Precompute UMAP dataframe
    umap_df = pd.DataFrame({
        "UMAP1": adata.obsm["X_umap"][:, 0],
        "UMAP2": adata.obsm["X_umap"][:, 1],
    })
    for col in cluster_cols:
        umap_df[col] = adata.obs[col].astype(str)

    # Default plot (first cluster resolution)
    fig_umap = px.scatter(
        umap_df,
        x="UMAP1", y="UMAP2",
        color=default_cluster,
        title=f"UMAP - {default_cluster}",
        opacity=0.7,
    )

    # Default bar plot (cells per cluster)
    cluster_counts = adata.obs[default_cluster].value_counts().reset_index()
    cluster_counts.columns = ["cluster", "n_cells"]
    fig_bar = px.bar(cluster_counts, x="cluster", y="n_cells", title="Cells per cluster")

    layout = html.Div([
        html.H3("Clustering and UMAP"),
        html.Label("Select clustering resolution:"),
        dcc.Dropdown(
            id="cluster-dropdown",
            options=[{"label": col, "value": col} for col in cluster_cols],
            value=default_cluster,
            clearable=False,
        ),
        dcc.Graph(id="umap-plot", figure=fig_umap),
        dcc.Graph(id="cluster-bar", figure=fig_bar)
    ])

    return layout