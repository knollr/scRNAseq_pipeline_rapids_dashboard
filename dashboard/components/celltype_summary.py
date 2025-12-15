"""
components/celltype_summary.py
----------------------------------
Cell type annotation summary (CellTypist output).

UMAP + stacked cell-type composition per source_file | Sample_Name
Supports AnnData and MuData input via data object.
"""

import numpy as np
import pandas as pd
import plotly.express as px
from dash import html, dcc, Input, Output


def _detect_models(adata):
    """Return list of CellTypist model prefixes."""
    prefixes = []
    for c in adata.obs.columns:
        if c.endswith("_predicted_labels"):
            prefixes.append(c.replace("_predicted_labels", ""))
    return sorted(prefixes)


def layout(data):
    """Return layout for CellTypist summaries."""
    adata = data.adata
    models = _detect_models(adata)
    if not models:
        models = ["No CellTypist models found"]

    dropdown_style = {"width": "100%", "marginBottom": "1rem"}

    return html.Div(
        [
            html.H3("🧠 Cell Type Annotation Summary", style={"marginBottom": "1rem"}),

            html.Div(
                [
                    html.Label("CellTypist model:", style={"font-weight": "bold"}),
                    dcc.Dropdown(
                        id="celltype-model-dropdown",
                        options=[{"label": m, "value": m} for m in models],
                        value=models[0],
                        clearable=False,
                        style=dropdown_style,
                    ),
                ]
            ),

            # --- Top row: UMAP (QC-style width) ---
            html.Div(
                [
                    dcc.Graph(id="celltype-umap"),
                    html.Div(),  # empty placeholder to keep 2-column layout
                ],
                style={
                    "display": "grid",
                    "gridTemplateColumns": "1fr 1fr",
                    "gap": "1.5rem",
                    "marginBottom": "1.5rem",
                },
            ),

            # --- Bottom row: wide stacked bar ---
            html.Div(
                [dcc.Graph(id="celltype-barplot")],
                style={
                    "display": "grid",
                    "gridTemplateColumns": "1fr",
                },
            ),
        ],
        style={"margin": "2rem"},
    )


def register_callbacks(app, data):
    adata = data.adata

    @app.callback(
        Output("celltype-umap", "figure"),
        Output("celltype-barplot", "figure"),
        Input("celltype-model-dropdown", "value"),
    )
    def update_celltype_plots(model_prefix):

        label_col = f"{model_prefix}_predicted_labels"
        if label_col not in adata.obs.columns:
            return (
                px.scatter(title="No cell type labels found"),
                px.bar(title="No data"),
            )

        df = adata.obs.copy()

        # --------------------------------------------------
        # Detect UMAP (MuData preferred)
        # --------------------------------------------------
        if hasattr(data, "mdata") and data.mdata is not None:
            if "X_wnn" in data.mdata.obsm:
                umap = data.mdata.obsm["X_wnn"][:, :2]
                idx = data.mdata.obs.index
            elif "X_umap" in data.mdata.obsm:
                umap = data.mdata.obsm["X_umap"][:, :2]
                idx = data.mdata.obs.index
            else:
                umap = None
        else:
            if "X_umap" in adata.obsm:
                umap = adata.obsm["X_umap"][:, :2]
                idx = adata.obs.index
            else:
                umap = None

        if umap is None:
            fig_umap = px.scatter(title="No UMAP found")
        else:
            umap_df = pd.DataFrame(umap, columns=["UMAP1", "UMAP2"], index=idx)

            # sampling
            n = len(umap_df)
            if n > 15000:
                sel = np.random.default_rng(0).choice(n, 15000, replace=False)
                plot_df = umap_df.iloc[sel]
            else:
                plot_df = umap_df

            plot_df[label_col] = df.loc[plot_df.index, label_col].astype(str)

            # consistent color map
            labels = plot_df[label_col].unique()
            palette = px.colors.qualitative.Plotly
            color_map = dict(zip(labels, (palette * 10)[: len(labels)]))

            fig_umap = px.scatter(
                plot_df,
                x="UMAP1",
                y="UMAP2",
                color=label_col,
                color_discrete_map=color_map,
                opacity=0.7,
                render_mode="webgl",
                template="plotly_white",
                title=f"UMAP colored by CellTypist: {model_prefix}",
            )
            fig_umap.update_traces(marker={"size": 3})
            fig_umap.update_layout(
                height=500,
                margin=dict(l=40, r=40, t=80, b=40),
            )

        # --------------------------------------------------
        # Stacked bar: source_file | Sample_Name
        # --------------------------------------------------
        if {"source_file", "Sample_Name"}.issubset(df.columns):
            df["source_sample"] = (
                df["source_file"].astype(str)
                + " | "
                + df["Sample_Name"].astype(str)
            )

            bar_df = (
                df.groupby(["source_sample", label_col])
                .size()
                .reset_index(name="Cell count")
            )

            fig_bar = px.bar(
                bar_df,
                x="source_sample",
                y="Cell count",
                color=label_col,
                color_discrete_map=color_map,
                barmode="stack",
                template="plotly_white",
                title="Cell type composition per source_file | Sample_Name",
            )

            fig_bar.update_layout(
                barnorm="fraction",
                height=450,
                margin=dict(l=40, r=40, t=80, b=120),
                xaxis_title="Source file | Sample name",
                yaxis=dict(
                    title="Fraction of cells",
                    tickformat=".0%",
                ),
            )
        else:
            fig_bar = px.bar(
                title="source_file or Sample_Name missing",
                template="plotly_white",
            )

        return fig_umap, fig_bar
