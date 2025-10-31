"""
components/qc_plots.py
----------------------------------
QC Overview dashboard component

Plots:
1️⃣ Violin plot for QC metrics per source_file (n_counts, n_genes, pct_counts_mt)
2️⃣ Bar plot: number of cells per source_file
3️⃣ UMAP: colored by source_file or source_file+Sample_Name
4️⃣ Stacked bar plot: cells per source_file, stacked by Sample_Name
"""

import numpy as np
import pandas as pd
import plotly.express as px
from dash import html, dcc, Input, Output


def layout(adata):
    """Return layout for QC visualization with dropdowns."""
    # detect metrics
    qc_metrics = [x for x in ["n_counts", "n_genes", "pct_counts_mt"] if x in adata.obs.columns]
    if not qc_metrics:
        qc_metrics = ["n_genes"]

    # responsive dropdown width
    dropdown_style = {"width": "60%", "marginBottom": "1rem"}

    return html.Div(
        [
            html.H3("📊 QC Overview", style={"marginBottom": "1rem"}),

            html.Div(
                [
                    html.Label("Select QC metric for violin plot:"),
                    dcc.Dropdown(
                        id="qc-metric-dropdown",
                        options=[{"label": m, "value": m} for m in qc_metrics],
                        value=qc_metrics[0],
                        clearable=False,
                        style=dropdown_style,
                    ),
                ]
            ),

            html.Div(
                [
                    html.Label("UMAP coloring:"),
                    dcc.Dropdown(
                        id="umap-color-dropdown",
                        options=[
                            {"label": "source_file", "value": "source_file"},
                            {"label": "source_file + Sample_Name", "value": "source_sample"},
                        ],
                        value="source_file",
                        clearable=False,
                        style={"width": "70%", "marginBottom": "1rem"},
                    ),
                ]
            ),

            html.Div(
                [
                    dcc.Graph(id="qc-violin"),
                    dcc.Graph(id="qc-bar-source"),
                    dcc.Graph(id="qc-umap"),
                    dcc.Graph(id="qc-bar-sample"),
                ],
                style={
                    "display": "grid",
                    "gridTemplateColumns": "1fr 1fr",
                    "gridTemplateRows": "auto auto",
                    "gap": "1.5rem",
                    "justifyItems": "center",
                },
            ),
        ],
        style={"margin": "2rem"},
    )


def register_callbacks(app, adata):
    """Register Dash callbacks for QC visualizations."""

    @app.callback(
        Output("qc-violin", "figure"),
        Output("qc-bar-source", "figure"),
        Output("qc-umap", "figure"),
        Output("qc-bar-sample", "figure"),
        Input("qc-metric-dropdown", "value"),
        Input("umap-color-dropdown", "value"),
    )
    def update_qc_plots(qc_metric, umap_choice):
        df = adata.obs.copy()

        # ===== Violin Plot =====
        if "source_file" in df.columns and qc_metric in df.columns:
            fig_violin = px.violin(
                df,
                x="source_file",
                y=qc_metric,
                box=True,
                points=False,
                title=f"Distribution of {qc_metric} per source_file",
                template="plotly_white",
            )
            fig_violin.update_traces(width=0.4)
            fig_violin.update_layout(
                xaxis_title="Source file",
                yaxis_title=qc_metric,
                margin=dict(l=40, r=40, t=80, b=40),
                height=400,
            )
        else:
            fig_violin = px.scatter(title="QC metric or source_file not found")

        # ===== Bar Plot: Cells per Source =====
        if "source_file" in df.columns:
            cell_counts = df["source_file"].value_counts().sort_index().reset_index()
            cell_counts.columns = ["source_file", "Cell count"]

            fig_bar_source = px.bar(
                cell_counts,
                x="source_file",
                y="Cell count",
                title="Cells per source_file",
                text_auto=True,
                template="plotly_white",
            )
            fig_bar_source.update_layout(height=400, margin=dict(l=40, r=40, t=80, b=40))
        else:
            fig_bar_source = px.bar(title="source_file column missing")

        # ===== UMAP =====
        if "X_umap" not in adata.obsm.keys():
            fig_umap = px.scatter(title="No UMAP found")
        else:
            umap_arr = adata.obsm["X_umap"][:, :2]
            umap_df = pd.DataFrame(umap_arr, columns=["UMAP1", "UMAP2"], index=adata.obs.index)

            n_total = umap_df.shape[0]
            umap_sample_limit = 15000
            if n_total > umap_sample_limit:
                sampled_pos = np.random.default_rng(0).choice(n_total, size=umap_sample_limit, replace=False)
                sampled_labels = umap_df.index[sampled_pos]
                umap_plot_df = umap_df.loc[sampled_labels]
            else:
                umap_plot_df = umap_df

            # --- Color logic ---
            if umap_choice == "source_sample" and {"source_file", "Sample_Name"}.issubset(df.columns):
                umap_plot_df["color"] = (
                    df.loc[umap_plot_df.index, "source_file"].astype(str)
                    + " | "
                    + df.loc[umap_plot_df.index, "Sample_Name"].astype(str)
                )
                title = "UMAP colored by source_file + Sample_Name"
            else:
                if "source_file" in df.columns:
                    umap_plot_df["color"] = df.loc[umap_plot_df.index, "source_file"].astype(str)
                    title = "UMAP colored by source_file"
                else:
                    umap_plot_df["color"] = "unknown"
                    title = "UMAP (no source_file found)"

            fig_umap = px.scatter(
                umap_plot_df,
                x="UMAP1",
                y="UMAP2",
                color="color",
                title=f"{title} (showing {umap_plot_df.shape[0]} / {n_total} cells)",
                opacity=0.7,
                render_mode="webgl",
                template="plotly_white",
            )
            fig_umap.update_traces(marker=dict(size=3, line=dict(width=0)))
            fig_umap.update_layout(height=500, margin=dict(l=40, r=40, t=80, b=40))

        # ===== Stacked Bar: Cells per Sample (per source) =====
        if {"source_file", "Sample_Name"}.issubset(df.columns):
            stacked = df.groupby(["source_file", "Sample_Name"]).size().reset_index(name="Cell count")
            fig_bar_sample = px.bar(
                stacked,
                x="source_file",
                y="Cell count",
                color="Sample_Name",
                title="Cells per source_file (stacked by Sample_Name)",
                template="plotly_white",
            )
            fig_bar_sample.update_layout(barmode="stack", height=400, margin=dict(l=40, r=40, t=80, b=40))
        else:
            fig_bar_sample = px.bar(title="Sample_Name or source_file missing")

        return fig_violin, fig_bar_source, fig_umap, fig_bar_sample
