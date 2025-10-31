"""
components/celltype_summary.py
----------------------------------
Cell type annotation summary (CellTypist output).
"""

import plotly.express as px
from dash import html, dcc, Input, Output


def _detect_models(adata):
    """Return list of CellTypist model prefixes."""
    obs_cols = adata.obs.columns
    prefixes = set()
    for c in obs_cols:
        if c.endswith("_predicted_labels"):
            prefixes.add(c.replace("_predicted_labels", ""))
    return sorted(prefixes)


def layout(adata):
    """Return layout for CellTypist summaries."""
    models = _detect_models(adata)
    if not models:
        models = ["No CellTypist models found"]

    return html.Div(
        [
            html.H3("🧠 Cell Type Annotation Summary"),
            dcc.Dropdown(
                id="model-dropdown",
                options=[{"label": m, "value": m} for m in models],
                value=models[0],
                clearable=False,
                style={"width": "50%"},
            ),
            html.Div(
                [
                    dcc.Graph(id="celltype-umap"),
                    dcc.Graph(id="conf-violin"),
                ],
                style={"display": "grid", "gridTemplateColumns": "1fr 1fr", "gap": "1rem"},
            ),
        ]
    )


def register_callbacks(app, adata):
    """Register callbacks for CellTypist annotation plots."""

    @app.callback(
        [Output("celltype-umap", "figure"), Output("conf-violin", "figure")],
        Input("model-dropdown", "value"),
    )
    def update_celltype_plots(model_prefix):
        df = adata.obs.copy()

        if f"{model_prefix}_predicted_labels" not in df.columns:
            return px.scatter(title="No celltype labels found"), px.violin(title="No data")

        df["UMAP1"] = adata.obsm["X_umap"][:, 0]
        df["UMAP2"] = adata.obsm["X_umap"][:, 1]

        label_col = f"{model_prefix}_predicted_labels"
        conf_col = f"{model_prefix}_conf_score"

        # --- UMAP by predicted labels ---
        fig_umap = px.scatter(
            df,
            x="UMAP1",
            y="UMAP2",
            color=label_col,
            title=f"UMAP colored by CellTypist model: {model_prefix}",
            opacity=0.7,
            render_mode="webgl",
        )
        fig_umap.update_traces(marker=dict(size=3, line=dict(width=0)))

        # --- Violin of confidence scores ---
        if conf_col in df.columns:
            fig_conf = px.violin(
                df,
                y=conf_col,
                x=label_col,
                color=label_col,
                title=f"Confidence score distribution per cell type ({model_prefix})",
                box=True,
                points=False,
            )
        else:
            fig_conf = px.violin(title="Confidence scores not available")

        return fig_umap, fig_conf
