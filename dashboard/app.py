# dashboard/app.py

import dash
from dash import Input, Output
import scanpy as sc
from components.clustering_umap import clustering_umap_component
import h5py

app = dash.Dash(__name__)

adata = sc.read_h5ad("/pipeline/results/run_2025-10-28_15-52-06/merged_annotated.h5ad")

app.layout = clustering_umap_component(adata)

@app.callback(
    [Output("umap-plot", "figure"),
     Output("cluster-bar", "figure")],
    Input("cluster-dropdown", "value")
)
def update_cluster_plots(cluster_col):
    import plotly.express as px
    import pandas as pd

    umap_df = pd.DataFrame({
        "UMAP1": adata.obsm["X_umap"][:, 0],
        "UMAP2": adata.obsm["X_umap"][:, 1],
        "cluster": adata.obs[cluster_col].astype(str)
    })

    fig_umap = px.scatter(
        umap_df,
        x="UMAP1", y="UMAP2",
        color="cluster",
        title=f"UMAP - {cluster_col}",
        opacity=0.7,
    )

    cluster_counts = adata.obs[cluster_col].value_counts().reset_index()
    cluster_counts.columns = ["cluster", "n_cells"]
    fig_bar = px.bar(cluster_counts, x="cluster", y="n_cells", title="Cells per cluster")

    return fig_umap, fig_bar

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=8050)
