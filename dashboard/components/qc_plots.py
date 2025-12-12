"""
components/qc_plots.py
----------------------------------
QC Overview dashboard component

2x2 QC plots:

A: Violin QC metrics
B: Cells per source_file OR cells per source_file stacked by Sample_Name
C: UMAP (color by source_file, source_file+Sample_Name, or QC metric)
D: Doublet QC panel (from pre-filtered dataframe)
"""

import numpy as np
import pandas as pd
import plotly.express as px
from dash import html, dcc, Input, Output


def layout(data):
    """Return layout for QC visualization with dropdowns."""
    adata = data.adata  # always the RNA AnnData (works for both modes)
    qc_metrics = [
        x
        for x in [
            "total_counts",
            "log1p_total_counts",
            "n_genes",
            "log1p_n_genes_by_counts",
            "pct_counts_mt",
            "pct_counts_ribo",
            "pct_counts_hb",
        ]
        if x in adata.obs.columns
    ]
    if not qc_metrics:
        qc_metrics = ["n_genes"]

    dropdown_style = {"width": "100%", "marginBottom": "1rem"}

    return html.Div(
        [
            html.H3("📊 QC Overview", style={"marginBottom": "1rem"}),

            # -------------------------------
            # DROPDOWNS (each on its own row)
            # -------------------------------
            html.Div(
                [
                    html.Label("1️⃣ - Violin QC metric:", style={"font-weight": "bold"}),
                    dcc.Dropdown(
                        id="qc-metric-dropdown",
                        options=[{"label": m, "value": m} for m in qc_metrics],
                        value=qc_metrics[0],
                        clearable=False,
                        style=dropdown_style,
                    ),
                ],
                style={"marginBottom": "10px"},
            ),

            html.Div(
                [
                    html.Label("2️⃣ – Cell Count View:", style={"font-weight": "bold"}),
                    dcc.Dropdown(
                        id="qc-bar-mode-dropdown",
                        options=[
                            {"label": "Cells per source_file", "value": "simple"},
                            {
                                "label": "Cells per source_file (stacked by Sample_Name)",
                                "value": "stacked",
                            },
                        ],
                        value="simple",
                        clearable=False,
                        style=dropdown_style,
                    ),
                ],
                style={"marginBottom": "10px"},
            ),

            html.Div(
                [
                    html.Label("3️⃣ – UMAP Coloring:", style={"font-weight": "bold"}),
                    dcc.Dropdown(
                        id="umap-color-dropdown",
                        options=[
                            {"label": "source_file", "value": "source_file"},
                            {
                                "label": "source_file + Sample_Name",
                                "value": "source_sample",
                            },
                        ]
                        + [{"label": m, "value": m} for m in qc_metrics],
                        value="source_file",
                        clearable=False,
                        style=dropdown_style,
                    ),
                ],
                style={"marginBottom": "20px"},
            ),
            
            html.Div(
                [
                    html.Label("4️⃣ – Doublet QC metric:", style={"font-weight": "bold"}),
                    dcc.Dropdown(
                        id="doublet-qc-metric-dropdown",
                        options=[
                            {"label": m, "value": m}
                            for m in [
                                "total_counts",
                                "log1p_total_counts",
                                "n_genes",
                                "log1p_n_genes_by_counts",
                                "pct_counts_mt",
                                "pct_counts_ribo",
                                "pct_counts_hb",
                                "doublet_score",
                            ]
                        ],
                        value="n_genes",
                        clearable=False,
                        style=dropdown_style,
                    ),
                ],
                style={"marginBottom": "20px"},
            ),


            # -------------------------------
            # 2×2 plot grid
            # -------------------------------
            html.Div(
                [
                    dcc.Graph(id="qc-violin"),
                    dcc.Graph(id="qc-bar-combined"),
                    dcc.Graph(id="qc-umap"),
                    dcc.Graph(id="qc-doublet-qc"),
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


def register_callbacks(app, data):
    adata = data.adata
    """Register Dash callbacks for QC visualizations."""

    @app.callback(
        Output("qc-violin", "figure"),
        Output("qc-bar-combined", "figure"),
        Output("qc-umap", "figure"),
        Output("qc-doublet-qc", "figure"),
        Input("qc-metric-dropdown", "value"),
        Input("qc-bar-mode-dropdown", "value"),
        Input("umap-color-dropdown", "value"),
        Input("doublet-qc-metric-dropdown", "value"),   # <-- ADD THIS

    )
    def update_qc_plots(qc_metric, bar_mode, umap_choice, doublet_qc_metric):

        df = adata.obs.copy()

        # consistent color map
        if "source_file" in df.columns:
            source_files = df["source_file"].unique()
            palette = px.colors.qualitative.Plotly
            repeats = int(np.ceil(len(source_files) / len(palette)))
            color_map = dict(zip(source_files, (palette * repeats)[: len(source_files)]))
        else:
            color_map = {}

        # ------------------------------------------------------------------
        # 1. VIOLIN PLOT
        # ------------------------------------------------------------------
        if "source_file" in df.columns and qc_metric in df.columns:
            fig_violin = px.violin(
                df,
                x="source_file",
                y=qc_metric,
                box=True,
                points=False,
                color="source_file",
                color_discrete_map=color_map,
                title=f"1️⃣ Distribution of {qc_metric} per source_file",
                template="plotly_white",
            )
            fig_violin.update_traces(width=0.8)
            fig_violin.update_layout(
                xaxis_title="Source file",
                yaxis_title=qc_metric,
                height=400,
                margin=dict(l=40, r=40, t=80, b=40),
            )
        else:
            fig_violin = px.scatter(title="source_file or metric missing")

        # ------------------------------------------------------------------
        # 2. COMBINED BAR PLOT
        # ------------------------------------------------------------------
        if bar_mode == "simple":
            if "source_file" in df.columns:
                cell_counts = (
                    df["source_file"].value_counts().sort_index().reset_index()
                )
                cell_counts.columns = ["source_file", "Cell count"]

                fig_bar = px.bar(
                    cell_counts,
                    x="source_file",
                    y="Cell count",
                    color="source_file",
                    color_discrete_map=color_map,
                    text_auto=True,
                    title="2️⃣ Cells per source_file",
                    template="plotly_white",
                )
                fig_bar.update_layout(
                    height=400, margin=dict(l=40, r=40, t=80, b=40)
                )
            else:
                fig_bar = px.bar(title="source_file missing")
        else:
            if {"source_file", "Sample_Name"}.issubset(df.columns):
                stacked = (
                    df.groupby(["source_file", "Sample_Name"])
                    .size()
                    .reset_index(name="Cell count")
                )

                fig_bar = px.bar(
                    stacked,
                    x="source_file",
                    y="Cell count",
                    color="Sample_Name",
                    title="2️⃣ Cells per source_file (stacked by Sample_Name)",
                    template="plotly_white",
                )
                fig_bar.update_layout(
                    barmode="stack",
                    height=400,
                    margin=dict(l=40, r=40, t=80, b=40),
                )
            else:
                fig_bar = px.bar(title="source_file or Sample_Name missing")

        # ------------------------------------------------------------------
        # 3. UMAP
        # ------------------------------------------------------------------
        # Detect embedding from AnnData or MuData
        if hasattr(data, "mdata") and data.mdata is not None:
            # --- MuData mode ---
            if "X_wnn" in data.mdata.obsm:
                umap_arr = data.mdata.obsm["X_wnn"][:, :2]
                embedding_index = data.mdata.obs.index
            elif "X_umap" in data.mdata.obsm:
                umap_arr = data.mdata.obsm["X_umap"][:, :2]
                embedding_index = data.mdata.obs.index
            else:
                umap_arr = None
        else:
            # --- AnnData mode ---
            if "X_umap" in adata.obsm:
                umap_arr = adata.obsm["X_umap"][:, :2]
                embedding_index = adata.obs.index
            else:
                umap_arr = None

        # If nothing found: STOP here
        if umap_arr is None:
            fig_umap = px.scatter(title="No UMAP / WNN UMAP found")
        else:
            umap_df = pd.DataFrame(
                umap_arr, columns=["UMAP1", "UMAP2"], index=embedding_index
            )

            # sampling
            n = len(umap_df)
            if n > 15000:
                idx = np.random.default_rng(0).choice(n, 15000, replace=False)
                plot_df = umap_df.iloc[idx].copy()
            else:
                plot_df = umap_df.copy()

            # ---- categorical coloring ----
            if umap_choice in ["source_file", "source_sample"]:
                if umap_choice == "source_sample" and {
                    "source_file",
                    "Sample_Name",
                }.issubset(df.columns):
                    plot_df["color"] = (
                        df.loc[plot_df.index, "source_file"].astype(str)
                        + " | "
                        + df.loc[plot_df.index, "Sample_Name"].astype(str)
                    )
                    map_colors = None
                    title = "UMAP colored by source_file + Sample_Name"
                else:
                    plot_df["color"] = df.loc[plot_df.index, "source_file"].astype(str)
                    map_colors = color_map
                    title = "3️⃣ UMAP colored by source_file"

                fig_umap = px.scatter(
                    plot_df,
                    x="UMAP1",
                    y="UMAP2",
                    color="color",
                    color_discrete_map=map_colors,
                    opacity=0.7,
                    template="plotly_white",
                    render_mode="webgl",
                    title=title,
                )
            else:
                # ---- numeric QC metric ----
                metric = umap_choice
                if metric not in df.columns:
                    fig_umap = px.scatter(title=f"Metric '{metric}' missing")
                else:
                    plot_df[metric] = df.loc[plot_df.index, metric]

                    fig_umap = px.scatter(
                        plot_df,
                        x="UMAP1",
                        y="UMAP2",
                        color=metric,
                        color_continuous_scale="Viridis",
                        opacity=0.8,
                        template="plotly_white",
                        render_mode="webgl",
                        title=f"3️⃣ UMAP colored by {metric}",
                    )

            fig_umap.update_traces(marker={"size": 3})
            fig_umap.update_layout(
                height=500, margin=dict(l=40, r=40, t=80, b=40)
            )

        # ------------------------------------------------------------------
        # 4. DOUBLETS PANEL (Plot 4)
        # ------------------------------------------------------------------
        if hasattr(data, "doublet_df") and data.doublet_df is not None:
            doublet_df = data.doublet_df.copy()

            # Ensure boolean
            if doublet_df["predicted_doublet"].dtype != bool:
                doublet_df["predicted_doublet"] = doublet_df["predicted_doublet"].astype(bool)

            # ----- QC metric dropdown options -----
            qc_violin_metrics = [
                "total_counts",
                "log1p_total_counts",
                "n_genes",
                "log1p_n_genes_by_counts",
                "pct_counts_mt",
                "pct_counts_ribo",
                "pct_counts_hb",
                "doublet_score",
            ]

            qc_metric_dd = (
                doublet_qc_metric if doublet_qc_metric in qc_violin_metrics else "n_genes"
            )

            # ----- Extract method & rate for subtitle -----
            dbl_method = (
                doublet_df["doublet_method"].iloc[0]
                if "doublet_method" in doublet_df.columns
                else "unknown"
            )

            dbl_rate = (
                doublet_df["doublet_rate"].iloc[0]
                if "doublet_rate" in doublet_df.columns
                else "n/a"
            )

            subtitle = f"Doublet detection with {dbl_method} (expected rate: {dbl_rate})"

            # ----- Violin plot -----
            fig_placeholder = px.violin(
                doublet_df,
                x="source_file",
                y=qc_metric_dd,
                color="predicted_doublet",
                box=True,
                points=False,
                template="plotly_white",
                title=f"4️⃣ {qc_metric_dd} by predicted doublet<br><sup>{subtitle}</sup>",
            )

            # Side-by-side violins
            fig_placeholder.update_layout(violinmode="group")
            fig_placeholder.update_traces(width=0.8)

            fig_placeholder.update_layout(
                height=450,
                margin=dict(l=40, r=40, t=80, b=40),
            )

        else:
            fig_placeholder = px.scatter(
                title="No doublet dataframe found",
                template="plotly_white",
                height=400,
            )

        return fig_violin, fig_bar, fig_umap, fig_placeholder
