#!/usr/bin/env python3
"""
export_cluster_markers.py

Compute marker genes for all clustering columns (Leiden/Louvain, all resolutions)
and export as a tidy CSV suitable for dashboards.

Handles both single-modal (AnnData, .h5ad) and multimodal (MuData, .h5mu) datasets.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import muon as mu

# ----------------------------
# Helper function: clustering from mdata to rna
# ----------------------------

def propagate_clustering_to_rna(mdata, rna_key="rna"):
    """
    Copy clustering columns (leiden_*, louvain_*) from mdata.obs
    into mdata.mod[rna_key].obs.

    This is required because clustering lives on MuData, but
    marker genes are computed on RNA AnnData.
    """
    if rna_key not in mdata.mod:
        raise ValueError(f"MuData has no modality '{rna_key}'")

    rna = mdata.mod[rna_key]

    cluster_cols = [
        c for c in mdata.obs.columns
        if c.startswith("leiden_") or c.startswith("louvain_")
    ]

    if not cluster_cols:
        raise ValueError("No clustering columns found in mdata.obs")

    print(f"🔁 Propagating {len(cluster_cols)} clustering columns to RNA modality")

    for col in cluster_cols:
        rna.obs[col] = mdata.obs.loc[rna.obs_names, col].values

    return rna


# ----------------------------
# Helper function: export markers
# ----------------------------
def export_all_cluster_markers(adata, output_file, top_n=50):
    """
    Compute marker genes for all clustering columns in .obs and export as tidy CSV.

    Detects all columns starting with 'leiden_' or 'louvain_' automatically.
    """
    cluster_cols = [c for c in adata.obs.columns if c.startswith("leiden_") or c.startswith("louvain_")]
    if not cluster_cols:
        raise ValueError("No clustering columns found in .obs (prefix 'leiden_' or 'louvain_').")

    all_rows = []

    for col in cluster_cols:
        cluster_type = "leiden" if col.startswith("leiden") else "louvain"
        resolution = float(col.split("_")[1]) if "_" in col else np.nan
        print(f"🧪 Computing markers for {col} ({cluster_type}, resolution={resolution})")

        sc.tl.rank_genes_groups(adata, groupby=col, method="wilcoxon", pts=True)

        groups = adata.obs[col].cat.categories if hasattr(adata.obs[col], 'cat') else np.unique(adata.obs[col])
        
        for grp in groups:
            for i in range(top_n):
                gene = adata.uns['rank_genes_groups']['names'][grp][i]
                logfc = adata.uns['rank_genes_groups']['logfoldchanges'][grp][i]
                pct_in = adata.uns['rank_genes_groups']['pts'][grp][i]
                pct_out = adata.uns['rank_genes_groups']['pts_rest'][grp][i]
                all_rows.append({
                    "cluster_col": col,
                    "cluster_type": cluster_type,
                    "resolution": resolution,
                    "cluster": grp,
                    "gene": gene,
                    "avg_logFC": logfc,
                    "pct_in_cluster": pct_in,
                    "pct_out_cluster": pct_out
                })

    df = pd.DataFrame(all_rows)
    df.to_csv(output_file, index=False)
    print(f"📁 Exported tidy cluster markers → {output_file} ({df.shape[0]} rows)")

# ----------------------------
# Main entry
# ----------------------------
def main(input_file, output_file, top_n=50):
    """
    Main entry point for Snakemake or CLI.
    Supports .h5ad (AnnData) and .h5mu (MuData) inputs.
    """
    if input_file.endswith(".h5mu"):
        mdata = mu.read_h5mu(input_file)
        # Compute markers only for RNA modality
        if "rna" not in mdata.mod:
            raise ValueError("MuData object has no RNA modality ('rna')")
        
        print("🧬 Preparing RNA modality for marker detection...")
        adata_rna = propagate_clustering_to_rna(mdata, rna_key="rna")
        
        print("🧬 Computing markers for RNA modality of MuData...")
        export_all_cluster_markers(adata_rna, output_file, top_n=top_n)

    elif input_file.endswith(".h5ad"):
        adata = sc.read(input_file)
        print("🧬 Computing markers for single-modal AnnData...")
        export_all_cluster_markers(adata, output_file, top_n=top_n)

    else:
        raise ValueError("Unsupported input format — must be .h5ad or .h5mu")

# ----------------------------
# CLI / Snakemake entry point
# ----------------------------
if __name__ == "__main__":
    try:
        input_file = snakemake.input[0]
        output_file = snakemake.output[0]
        top_n = getattr(snakemake.params, "top_n", 50)
    except NameError:
        # CLI execution
        parser = argparse.ArgumentParser(description="Export cluster marker genes to tidy CSV")
        parser.add_argument("--input", "-i", required=True, help="Input preprocessed AnnData (.h5ad) or MuData (.h5mu)")
        parser.add_argument("--output", "-o", required=True, help="Output CSV file")
        parser.add_argument("--top_n", type=int, default=50, help="Top N markers per cluster")
        args = parser.parse_args()
        input_file = args.input
        output_file = args.output
        top_n = args.top_n

    main(input_file, output_file, top_n)
