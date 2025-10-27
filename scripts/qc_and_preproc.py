#!/usr/bin/env python3
"""
qc_and_preproc.py

Preprocessing and QC script for single-cell datasets in a Snakemake pipeline.

Supports both:
- Single-modality (AnnData / .h5ad)
- Multimodal (MuData / .h5mu)

Features:
- RNA filtering, normalization, and dimensionality reduction
- CLR normalization for protein modality (CITE-seq)
- Weighted Nearest Neighbor (WNN) analysis for multimodal integration
- UMAP embedding and Leiden/Louvain clustering (CPU or GPU)
- Configurable storage of raw counts (adata.layers["counts"])
- Compatible with Snakemake and CLI
"""
#import scripts.silence_warnings
import argparse
import yaml
import scanpy as sc
import muon as mu
import numpy as np
import scipy.sparse as sp


# ----------------------------
# Helper: clustering
# ----------------------------
def run_clusterings(obj, resolutions, use_gpu=False, wnn_key=None, cfg=None):
    """
    Run Leiden and Louvain clustering for multiple resolutions.

    Automatically detects whether to use Scanpy (CPU) or RAPIDS (GPU),
    and supports clustering on multimodal WNN graphs via neighbors_key="wnn".

    Parameters
    ----------
    obj : AnnData or MuData
        Object to cluster
    resolutions : list of float
        List of clustering resolutions
    use_gpu : bool
        Use rapids-singlecell if True
    wnn_key : str or None
        Key of neighbors graph to use (e.g. "wnn" for multimodal)
    """
    is_mudata = hasattr(obj, "mod")
    
    clustering_cfg = cfg["clustering"]
    run_leiden = clustering_cfg.get("leiden", True)
    run_louvain = clustering_cfg.get("louvain", True)

    # Select backend
    if use_gpu:
        import rapids_singlecell as rsc
        if run_leiden:
            leiden_func = rsc.tl.leiden
        if run_louvain:
            louvain_func = rsc.tl.louvain
        backend = "⚡ RAPIDS (GPU)"
    else:
        if run_leiden:
            leiden_func = sc.tl.leiden
        if run_louvain:
            louvain_func = sc.tl.louvain
        backend = "🧠 Scanpy (CPU)"

    print(f"🏃 Running clustering with {backend}")
    if is_mudata and wnn_key:
        print(f"🧬 Using multimodal neighbors graph: '{wnn_key}'")

    for res in resolutions:
        if run_leiden:
            kwargs = {
                "resolution": res,
                "n_iterations": 2,
                "key_added": f"leiden_{res}" if not wnn_key else f"leiden_{res}_{wnn_key}",
                "neighbors_key": wnn_key
            }

            if not use_gpu:
                kwargs["flavor"] = "igraph"  # only for CPU, does not work with rsc
                kwargs["directed"] = False
            leiden_func(obj, **kwargs)

        if run_louvain:    
            louvain_func(
                obj,
                resolution=res,
                key_added=f"louvain_{res}" if not wnn_key else f"louvain_{res}_{wnn_key}",
                neighbors_key=wnn_key
            )


# ----------------------------
# RNA preprocessing
# ----------------------------
def filter_qc(adata, cfg):
    """
    Filter, normalize, and preprocess RNA data using Scanpy.

    Parameters
    ----------
    adata : AnnData
        RNA AnnData object
    cfg : dict
        Configuration dictionary

    Returns
    -------
    AnnData
        Preprocessed AnnData
    """
    p = cfg["params"]
    sc.pp.filter_cells(adata, min_genes=p["min_genes"])
    sc.pp.filter_genes(adata, min_cells=p["min_cells"])

    # Optionally store raw counts
    if cfg.get("store_raw_counts", False):
        print("💾 Storing raw counts in adata.layers['counts']")
        adata.layers["counts"] = adata.X.copy()

    # Normalization and feature selection
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
            
    print("🧹 Robust cleanup (CPU): handling sparse/dense/masked inputs...")

    # Case A: sparse matrix (Scipy)
    if sp.issparse(adata.X):
        # ensure float dtype in data array
        if not np.issubdtype(adata.X.data.dtype, np.floating):
            adata.X.data = adata.X.data.astype(np.float64)
        # Replace NaNs/Infs in the sparse data array
        if np.isnan(adata.X.data).any() or np.isinf(adata.X.data).any():
            print("⚠️ NaNs/Infs detected in sparse .X.data → replacing with 0")
            adata.X.data = np.nan_to_num(adata.X.data, nan=0.0, posinf=0.0, neginf=0.0)
        # Remove all-zero genes (columns)
        gene_counts = np.asarray(adata.X.sum(axis=0)).ravel()
        nonzero_genes = gene_counts > 0
        if not nonzero_genes.all():
            n_removed = (~nonzero_genes).sum()
            print(f"🧬 Removing {n_removed} all-zero genes (sparse)")
            adata = adata[:, nonzero_genes].copy()

    # Case B: dense / masked / array-like
    else:
        # coerce to numpy array (handles pandas-backed, masked arrays, etc.)
        X = np.asarray(adata.X)

        # If masked array, fill with 0
        if np.ma.isMaskedArray(X):
            X = X.filled(0.0)

        # Ensure floating dtype
        if not np.issubdtype(X.dtype, np.floating):
            X = X.astype(np.float64)
        else:
            X = X.astype(np.float64, copy=False)

        # Replace NaN / Inf if present
        if np.isnan(X).any() or np.isinf(X).any():
            print("⚠️ NaNs/Infs detected in dense .X → replacing with 0")
            X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

        # Remove all-zero genes (columns)
        gene_counts = (X > 0).sum(axis=0)
        nonzero_genes = gene_counts > 0
        if not nonzero_genes.all():
            n_removed = (~nonzero_genes).sum()
            print(f"🧬 Removing {n_removed} all-zero genes (dense)")
            adata = adata[:, nonzero_genes].copy()
        else:
            adata.X = X  # write cleaned dense matrix back  
            
    adata.layers["norm_data"] = adata.X.copy()
    
    sc.pp.highly_variable_genes(adata, n_top_genes=p["n_top_genes"])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=p["n_pcs"])

    # Neighborhood graph, UMAP, clustering
    sc.pp.neighbors(adata, n_neighbors=p["n_neighbors"], n_pcs=p["n_pcs"])
    sc.tl.umap(adata, n_components=p["umap_n_components"])
    run_clusterings(adata, resolutions=[0.2, 0.4, 0.6, 0.8, 1.0], cfg=cfg)
    return adata


# ----------------------------
# CPU and GPU pipelines
# ----------------------------
def run_cpu_pipeline(adata, cfg):
    print("🧠 Running CPU-based Scanpy preprocessing...")
    return filter_qc(adata, cfg)


def run_gpu_pipeline(adata, cfg):
    """
    GPU-accelerated preprocessing using rapids-singlecell.
    """
    import rapids_singlecell as rsc
    from cupyx.scipy.sparse import issparse
    import cupy as cp

    p = cfg["params"]
    print("⚡ Running GPU-based preprocessing with rapids-singlecell...")

    # Ensure .X is float32 (compatible with CuPy)
    if sp.issparse(adata.X):
        adata.X = adata.X.astype(np.float32)
    else:
        adata.X = np.array(adata.X, dtype=np.float32)

    rsc.get.anndata_to_GPU(adata)
    rsc.pp.filter_cells(adata, min_genes=p["min_genes"])
    rsc.pp.filter_genes(adata, min_cells=p["min_cells"])

    if cfg.get("store_raw_counts", False):
        print("💾 Storing raw counts in adata.layers['counts']")
        adata.layers["counts"] = adata.X.copy()
        rsc.get.anndata_to_CPU(adata, layers=["counts"])

    rsc.pp.normalize_total(adata, target_sum=1e4)
    rsc.pp.log1p(adata)
    
    # Remove NaN or invalid values directly on GPU
    print("🧹 Checking for NaN or empty features (GPU mode)...")
    # Detect and replace NaNs / infinities 
    if issparse(adata.X):
        # Only operate on the stored data
        if cp.isnan(adata.X.data).any():
            print("⚠️ NaNs detected → replacing with 0")
            adata.X.data = cp.nan_to_num(
                adata.X.data, nan=0.0, posinf=0.0, neginf=0.0
            )
    else:
        if cp.isnan(adata.X).any():
            print("⚠️ NaNs detected → replacing with 0")
            adata.X = cp.nan_to_num(
                adata.X, nan=0.0, posinf=0.0, neginf=0.0
            )
    # Remove all-zero genes (columns) 
    if issparse(adata.X):
        # Count non-zero entries per column
        X_pos = (adata.X > 0).astype(cp.float32)  # <-- important!
        nonzero_genes = cp.array(X_pos.sum(axis=0)).ravel() > 0
    else:
        nonzero_genes = cp.array((adata.X > 0).sum(axis=0)).ravel() > 0
    if not nonzero_genes.all():
        print(f"🧬 Removing {(~nonzero_genes).sum()} all-zero genes")
        # Indexing works with sparse or dense matrices
        adata = adata[:, cp.asnumpy(nonzero_genes)].copy()
        
    adata.layers["norm_data"] = adata.X.copy()    
    rsc.pp.highly_variable_genes(adata, n_top_genes=p["n_top_genes"])
    rsc.pp.scale(adata)
    rsc.tl.pca(adata, n_comps=p["n_pcs"])
    rsc.pp.neighbors(adata, n_neighbors=p["n_neighbors"], n_pcs=p["n_pcs"])
    rsc.tl.umap(adata)
    run_clusterings(adata, resolutions=[0.2, 0.4, 0.6, 0.8, 1.0], use_gpu=True, cfg=cfg)
    rsc.get.anndata_to_CPU(adata, convert_all = True)
    
    return adata


# ----------------------------
# Multimodal preprocessing
# ----------------------------
def preprocess_mdata(mdata, cfg, use_gpu):
    """
    Preprocess MuData object (RNA + protein multimodal).

    - RNA preprocessing
    - Protein CLR normalization
    - Weighted Nearest Neighbor (WNN) integration
    - UMAP and clustering on WNN representation
    """
    p = cfg["params"]

    # RNA modality
    print("🧬 Preprocessing RNA modality...")
    if use_gpu:
        mdata.mod["rna"] = run_gpu_pipeline(mdata.mod["rna"], cfg)
    else:
        mdata.mod["rna"] = run_cpu_pipeline(mdata.mod["rna"], cfg)

    remaining_cells = mdata.mod['rna'].obs_names
    mdata.mod['prot'] = mdata.mod['prot'][remaining_cells, :].copy()
    mdata.update()

    # Protein modality
    if "prot" in mdata.mod:
        from muon import prot as pt
        print("💧 Preprocessing protein modality (CLR normalization)...")
        prot = mdata.mod["prot"]
        prot.X = prot.X.astype("float64")
        pt.pp.clr(prot)
        sc.pp.scale(prot, max_value=10)
        sc.tl.pca(prot, mask_var=None, n_comps=p["n_pcs_prot"])
        sc.pp.neighbors(prot, n_neighbors=p["n_neighbors"], n_pcs=p["n_pcs_prot"])
        mdata.mod["prot"] = prot

    # Weighted Nearest Neighbor integration
    print("🔗 Computing WNN multimodal integration...")
    
    mu.pp.neighbors(
        mdata, 
        key_added="wnn",  # name of the WNN graph in mdata.obsp
        add_weights_to_modalities=True  # optional: add modality weights to each modality
    )   
    
    # UMAP and clustering on WNN
    print("🏃 Running UMAP + clustering on WNN graph...")
    if use_gpu:
        mu.tl.umap(mdata, neighbors_key="wnn", method="rapids")
    else:
        mu.tl.umap(mdata, neighbors_key="wnn")

    run_clusterings(mdata, resolutions=[0.2, 0.4, 0.6, 0.8, 1.0], use_gpu=use_gpu, wnn_key="wnn", cfg=cfg)
    return mdata


# ----------------------------
# Main entry
# ----------------------------
def main(input_file, output_file, cfg, use_gpu):

    if input_file.endswith(".h5mu"):
        mdata = mu.read_h5mu(input_file)
        mdata = preprocess_mdata(mdata, cfg, use_gpu)
        mdata.write_h5mu(output_file)
        print(f"✅ Multimodal preprocessing finished → {output_file}")

    elif input_file.endswith(".h5ad"):
        adata = sc.read(input_file)
        adata = run_gpu_pipeline(adata, cfg) if use_gpu else run_cpu_pipeline(adata, cfg)
        adata.write_h5ad(output_file)
        print(f"✅ Single-modal preprocessing finished → {output_file}")

    else:
        raise ValueError("Unsupported input format — must be .h5ad or .h5mu")


# ----------------------------
# CLI / Snakemake entry point
# ----------------------------
if __name__ == "__main__":
    try:
        # Snakemake execution
        input_file = snakemake.input[0]
        output_file = snakemake.output[0]
        cfg = snakemake.params.config  # already a dict
        use_gpu = snakemake.params.use_gpu
    except NameError:
        # CLI execution
        ap = argparse.ArgumentParser(description="QC + preprocessing for single-cell datasets")
        ap.add_argument("--input", "-i", required=True, help="Input .h5ad or .h5mu file")
        ap.add_argument("--output", "-o", required=True, help="Output file")
        ap.add_argument("--config", "-c", required=True, help="YAML configuration file")
        args = ap.parse_args()
        with open(args.config) as f:
            cfg = yaml.safe_load(f)
        input_file, output_file = args.input, args.output

    main(input_file, output_file, cfg, use_gpu)