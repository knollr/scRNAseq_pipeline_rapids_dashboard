"""
components/clustering_umap.py
----------------------------------
Clustering & expression dashboard.

2x2 layout:

P1: UMAP colored by clustering
P2: Cells per cluster (same colors as P1)
P3: Cluster markers (placeholder)
P4: Feature expression UMAP (RNA / Protein)
"""

import numpy as np
import pandas as pd
import plotly.express as px
from dash import html, dcc, Input, Output
import tempfile
import io
import base64
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scanpy as sc
import os
import plotly.express as px

# ---------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------
def get_top_markers(marker_df, cluster_col, top_n=3):
    """
    Extract top N marker genes per cluster
    using explicit cluster_col match.
    """

    df = marker_df.copy()

    df = df[df["cluster_col"] == cluster_col]

    if df.empty:
        raise ValueError(f"No markers found for {cluster_col}")

    # Ensure cluster is string for consistency
    df["cluster"] = df["cluster"].astype(str)

    # Top N genes per cluster by avg_logFC
    top = (
        df.sort_values("avg_logFC", ascending=False)
        .groupby("cluster", group_keys=False)
        .head(top_n)
    )

    genes_by_cluster = (
        top.groupby("cluster")["gene"]
        .apply(list)
        .to_dict()
    )

    return genes_by_cluster

def prepare_dotplot_df(adata_rna, genes_by_cluster, cluster_col):
    """
    Returns a dataframe suitable for Plotly dotplot:
    Columns: gene, cluster, mean_expr, fraction
    """
    import pandas as pd
    import numpy as np

    rows = []
    for cluster, genes in genes_by_cluster.items():
        cells_in_cluster = adata_rna.obs[cluster_col] == cluster
        expr_sub = adata_rna[cells_in_cluster, genes].X
        if not isinstance(expr_sub, np.ndarray):
            expr_sub = expr_sub.toarray()

        mean_expr = expr_sub.mean(axis=0)
        frac = (expr_sub > 0).sum(axis=0) / expr_sub.shape[0]

        for g, m, f in zip(genes, mean_expr, frac):
            rows.append({"gene": g, "cluster": str(cluster), "mean_expr": m, "fraction": f})

    df = pd.DataFrame(rows)
    return df




def plotly_dotplot(df, title="Dotplot"):
    """
    df: gene / cluster / mean_expr / fraction
    """
    fig = px.scatter(
        df,
        x="gene",
        y="cluster",
        size="fraction",
        color="mean_expr",
        color_continuous_scale="Viridis",
        size_max=12,
        text=None,
        template="plotly_white",
    )

    fig.update_layout(
        title=title,
        xaxis_tickangle=-90,
        xaxis=dict(tickfont=dict(size=10)),
        yaxis=dict(autorange="reversed"),
        height= max(500, 25 * df['cluster'].nunique()),
        margin=dict(l=80, r=20, t=50, b=120),
    )

    return fig

# ---------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------
def layout(data):
    adata = data.adata
    mdata = getattr(data, "mdata", None)

    if mdata is not None:
        obs = mdata.obs
    else:
        obs = adata.obs

    cluster_cols = [
        c for c in obs.columns
        if c.startswith("leiden") or c.startswith("louvain")
    ]
    if not cluster_cols:
        cluster_cols = ["no_cluster_found"]

    pbmc_marker_genes = [
        # T cells
        "CD3D", "CD3E","IL7R","LTB",
        # CD8 / cytotoxic
        "NKG7","GNLY","GZMB", "CD8A",
        # NK cells
        "KLRD1","FCGR3A", "KLRF1",
        # B cells
        "MS4A1","CD79A",
        # Monocytes (CD14+ / FCGR3A+)
        "LYZ","S100A8","S100A9","CTSD", "CD14", 
        # Dendritic cells
        "FCER1A","CST3",
        # Platelets
        "PPBP","PF4",
    ]

    has_protein = (
        hasattr(data, "mdata")
        and data.mdata is not None
        and "prot" in data.mdata.mod
    )

    dropdown_style = {"width": "100%", "marginBottom": "1rem"}

    controls = [
        html.Div(
            [
                html.Label("1️⃣+2️⃣ Clustering:", style={"font-weight": "bold"}),
                dcc.Dropdown(
                    id="cluster-dropdown",
                    options=[{"label": c, "value": c} for c in cluster_cols],
                    value=cluster_cols[0],
                    clearable=False,
                    style=dropdown_style,
                ),
            ]
        )
    ]

    if has_protein:
        controls.append(
            html.Div(
                [
                    html.Label("4️⃣ Modality:", style={"font-weight": "bold"}),
                    dcc.Dropdown(
                        id="modality-dropdown",
                        options=[
                            {"label": "RNA", "value": "rna"},
                            {"label": "Protein", "value": "protein"},
                        ],
                        value="rna",  # ✅ default RNA
                        clearable=False,
                        style=dropdown_style,
                    ),
                ]
            )
        )

    controls.extend(
        [
            html.Div(
                [
                    html.Label("4️⃣ Select feature (type search possible):", style={"font-weight": "bold"}),
                    dcc.Dropdown(
                        id="feature-dropdown",
                        options=[{"label": g, "value": g} for g in pbmc_marker_genes],
                        value=pbmc_marker_genes[0],
                        clearable=False,
                        searchable=True,  # ✅ search enabled
                        style=dropdown_style,
                    ),
                ]
            ),
            html.Div(
                [
                    dcc.Checklist(
                        id="log1p-toggle",
                        options=[{"label": " log1p transform", "value": "log1p"}],
                        value=[],
                    )
                ],
                style={"marginBottom": "1.5rem"},
            ),
        ]
    )

    return html.Div(
        [
            html.H3("🧬 Clustering & Expression Overview", style={"marginBottom": "1rem"}),
            html.Div(controls),
            html.Div(
                [
                    dcc.Graph(id="cluster-umap"),
                    dcc.Graph(id="cluster-barplot"),
                    dcc.Graph(id="cluster-marker-placeholder"),
                    dcc.Graph(id="feature-umap"),
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


# ---------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------
def register_callbacks(app, data):
    adata = data.adata

    # ------------------------------------------------------------------
    # Feature dropdown update
    # ------------------------------------------------------------------
    @app.callback(
        Output("feature-dropdown", "options"),
        Output("feature-dropdown", "value"),
        Input("modality-dropdown", "value"),
        prevent_initial_call=True,
    )
    def update_feature_dropdown(modality):

        pbmc_marker_genes = [
            # T cells
            "CD3D", "CD3E","IL7R","LTB",  
            # CD8 / cytotoxic
            "NKG7","GNLY","GZMB", "CD8A",
            # NK cells
            "KLRD1","FCGR3A", "KLRF1", 
            # B cells
            "MS4A1","CD79A",
            # Monocytes (CD14+ / FCGR3A+)
            "LYZ","S100A8","S100A9","CTSD", "CD14", 
            # Dendritic cells
            "FCER1A","CST3",
            # Platelets
            "PPBP","PF4",
        ]


        if modality == "protein":
            prot = data.mdata.mod["prot"]
            opts = [{"label": g, "value": g} for g in prot.var_names]
            return opts, prot.var_names[0]

        opts = [{"label": g, "value": g} for g in pbmc_marker_genes]
        return opts, pbmc_marker_genes[0]

    # ------------------------------------------------------------------
    # Main plots
    # ------------------------------------------------------------------
    @app.callback(
        Output("cluster-umap", "figure"),
        Output("cluster-barplot", "figure"),
        Output("cluster-marker-placeholder", "figure"),
        Output("feature-umap", "figure"),
        Input("cluster-dropdown", "value"),
        Input("modality-dropdown", "value"),
        Input("feature-dropdown", "value"),
        Input("log1p-toggle", "value"),
    )
    def update_plots(cluster_col, modality, feature, log1p_toggle):

        df = adata.obs.copy()
        if hasattr(data, "mdata") and data.mdata is not None:
            df = data.mdata.obs.copy()
        else:
            df = adata.obs.copy()
        log1p = "log1p" in log1p_toggle

        # ------------------------------------------------------------------
        # UMAP detection (QC-style)
        # ------------------------------------------------------------------
        if hasattr(data, "mdata") and data.mdata is not None:
            if "X_wnn" in data.mdata.obsm:
                emb = data.mdata.obsm["X_wnn"][:, :2]
                emb_index = data.mdata.obs.index
            elif "X_umap" in data.mdata.obsm:
                emb = data.mdata.obsm["X_umap"][:, :2]
                emb_index = data.mdata.obs.index
            else:
                emb = None
        else:
            if "X_umap" in adata.obsm:
                emb = adata.obsm["X_umap"][:, :2]
                emb_index = adata.obs.index
            else:
                emb = None

        if emb is None:
            empty = px.scatter(title="No UMAP / WNN found")
            return empty, empty, empty, empty

        umap_df = pd.DataFrame(emb, columns=["UMAP1", "UMAP2"], index=emb_index)

        # sampling
        n = len(umap_df)
        if n > 15000:
            idx = np.random.default_rng(0).choice(n, 15000, replace=False)
            plot_df = umap_df.iloc[idx].copy()
        else:
            plot_df = umap_df.copy()

        # ------------------------------------------------------------------
        # Shared cluster colors
        # ------------------------------------------------------------------
        if cluster_col in df.columns:
            clusters = df[cluster_col].astype(str).unique()
            palette = px.colors.qualitative.Plotly
            repeats = int(np.ceil(len(clusters) / len(palette)))
            cluster_colors = dict(
                zip(clusters, (palette * repeats)[: len(clusters)])
            )
        else:
            cluster_colors = {}

        # ------------------------------------------------------------------
        # P1: clustering UMAP
        # ------------------------------------------------------------------
        plot_df["cluster"] = df.loc[plot_df.index, cluster_col].astype(str)

        fig_cluster_umap = px.scatter(
            plot_df,
            x="UMAP1",
            y="UMAP2",
            color="cluster",
            color_discrete_map=cluster_colors,
            opacity=0.7,
            render_mode="webgl",
            template="plotly_white",
            title=f"1️⃣ UMAP colored by {cluster_col}",
        )
        fig_cluster_umap.update_traces(marker={"size": 3})
        fig_cluster_umap.update_layout(height=500)

        # ------------------------------------------------------------------
        # P2: cells per cluster (same colors)
        # ------------------------------------------------------------------
        counts = (
            df[cluster_col]
            .astype(str)
            .value_counts()
            .sort_values(ascending=False)
            .reset_index()
        )
        counts.columns = ["Cluster", "Cell count"]

        fig_bar = px.bar(
            counts,
            x="Cluster",
            y="Cell count",
            color="Cluster",
            color_discrete_map=cluster_colors,
            text_auto=True,
            template="plotly_white",
            title=f"2️⃣ Cells per {cluster_col}",
        )
        fig_bar.update_layout(height=500)

        # ------------------------------------------------------------------
        # P3: cluster marker dotplot (PRECOMPUTED markers)
        # ------------------------------------------------------------------
        try:
            # --- select RNA AnnData ---
            if hasattr(data, "mdata") and data.mdata is not None:
                adata_rna = data.mdata.mod["rna"]
                obs = data.mdata.obs
            else:
                adata_rna = data.adata
                obs = data.adata.obs

            # --- ensure clustering exists on RNA adata ---
            if cluster_col not in adata_rna.obs:
                adata_rna.obs[cluster_col] = obs[cluster_col].astype(str)

            # --- extract top markers ---
            genes_by_cluster = get_top_markers(
                marker_df=data.marker_df,
                cluster_col=cluster_col,
                top_n=3,
            )

            df_dot = prepare_dotplot_df(adata_rna, genes_by_cluster, cluster_col)
            fig_marker = plotly_dotplot(df_dot, title=f"3️⃣ Top markers per {cluster_col}")


        except Exception as e:
            fig_marker = px.scatter(
                title=f"3️⃣ Marker dotplot unavailable: {e}",
                template="plotly_white",
            )
            fig_marker.update_layout(height=500)

        # ------------------------------------------------------------------
        # P4: feature expression UMAP
        # ------------------------------------------------------------------
        if modality == "protein":
            prot = data.mdata.mod["prot"]
            expr = prot[plot_df.index, feature].X
        else:
            expr = adata[plot_df.index, feature].X

        if not isinstance(expr, np.ndarray):
            expr = expr.toarray()

        expr = expr.ravel()
        if log1p:
            expr = np.log1p(expr)

        plot_df["expr"] = expr

        fig_feat = px.scatter(
            plot_df,
            x="UMAP1",
            y="UMAP2",
            color="expr",
            color_continuous_scale="Plasma" if modality == "protein" else "Viridis",
            opacity=0.8,
            render_mode="webgl",
            template="plotly_white",
            title=f"4️⃣ {modality.upper()} expression: {feature}",
        )
        fig_feat.update_traces(marker={"size": 3})
        fig_feat.update_layout(height=500)

        return fig_cluster_umap, fig_bar, fig_marker, fig_feat
