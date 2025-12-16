#!/usr/bin/env python3
"""
celltypist_annotation.py

Annotate cell types using CellTypist in a Snakemake pipeline.

Supports:
- Single-modality (.h5ad)
- Multimodal (.h5mu, annotates RNA modality)
- Multiple models (from config)
- Optional GPU acceleration (RAPIDS backend)

Stores per-model:
- predicted labels
- confidence scores
- majority-voting labels (if enabled)
"""
import scripts.silence_warnings
import argparse
import os
import yaml
import scanpy as sc
import muon as mu
import celltypist
from celltypist import models
import numpy as np
import scipy.sparse as sp


# ----------------------------
# Helpers
# ----------------------------
def run_celltypist_annotation(adata, cfg, use_gpu=False):
    """
    Run CellTypist annotation for one or multiple models.

    Parameters
    ----------
    adata : AnnData
        RNA modality data for annotation
    cfg : dict
        CellTypist configuration
    use_gpu : bool
        Whether to use RAPIDS backend
    """
    models_list = cfg.get("models", [])
    if isinstance(models_list, str):
        models_list = [models_list]
    if not models_list:
        raise ValueError("No models specified under config['celltypist']['models'].")

    majority_voting = cfg.get("majority_voting", True)
    force_update = cfg.get("force_update", False)
    backend = "RAPIDS (GPU)" if use_gpu else "CPU"

    # Info about model path
    model_storage = models.models_path
    print(f"🧬 CellTypist model directory: {model_storage}")
    print(f"🧠/⚡ Backend: {backend} | 🗳️  Majority voting: {majority_voting}")

    for model_name in models_list:
        print(f"🖋️  Annotating with model: {model_name}")

        # Ensure model is present locally
        models.download_models(model=model_name, force_update=force_update)

        # Run annotation
        predictions = celltypist.annotate(
            adata,
            model=model_name,
            majority_voting=majority_voting,
            use_GPU = True if use_gpu else False
        )

        tag_prefix = f"{model_name}_"
        pred_adata = predictions.to_adata(prefix=tag_prefix)
        cols_to_add = [c for c in pred_adata.obs.columns if c.startswith(tag_prefix)]
        for col in cols_to_add:
            adata.obs[col] = pred_adata.obs[col]

    print(f"✅ Added CellTypist predictions to adata.obs: {cols_to_add}")

    adata.uns["celltypist_summary"] = {
        "models": models_list,
        "use_gpu": use_gpu,
        "majority_voting": majority_voting,
        "model_dir": model_storage,
        "version": celltypist.__version__,
    }

    return adata


def prepare_adata_for_celltypist(adata):
    """
    Prepare the expression matrix for CellTypist:
    - Prefer normalized layer if exists
    - Otherwise normalize counts and log1p
    - Ensure dtype float64
    - Remove NaN/Inf values
    """
    print("🧹 Preparing expression matrix for CellTypist...")

    if "norm_data" in adata.layers: # and use_gpu == False:
        print("✅ Using adata.layers['norm_data']")
        adata.X = adata.layers["norm_data"].copy().astype(np.float64)
    elif "counts" in adata.layers:
        print("🔄 Using adata.layers['counts'] and normalizing (fallback)")
        adata.X = adata.layers["counts"].copy().astype(np.float64)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("⚠️ No 'norm_data' or 'counts' layer found — using current .X as-is")
        adata.X = np.array(adata.X, dtype=np.float64)

    # Ensure float64 and replace NaN/Inf
    if sp.issparse(adata.X):
        if not np.issubdtype(adata.X.data.dtype, np.floating):
            adata.X.data = adata.X.data.astype(np.float64)
        if np.isnan(adata.X.data).any() or np.isinf(adata.X.data).any():
            print("⚠️ NaN/Inf detected in sparse .X.data → replacing with 0")
            adata.X.data = np.nan_to_num(adata.X.data, nan=0.0, posinf=0.0, neginf=0.0)
    else:
        X = np.asarray(adata.X)
        if np.ma.isMaskedArray(X):
            X = X.filled(0.0)
        if not np.issubdtype(X.dtype, np.floating):
            X = X.astype(np.float64)
        if np.isnan(X).any() or np.isinf(X).any():
            print("⚠️ NaN/Inf detected in dense .X → replacing with 0")
            X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
        adata.X = X

    totals = np.asarray(adata.X.sum(axis=1)).ravel()
    median_total = float(np.median(totals))
    min_x = float(np.min(adata.X))
    max_x = float(np.max(adata.X))

    print(f"ℹ️ Median total counts per cell (after log1p): {median_total:.2f}")
    print(f"ℹ️ X value range: min={min_x:.4f}, max={max_x:.4f}")

    if not np.isclose(median_total, 1e4, rtol=0.25):
        print(f"⚠️ Median total count deviates from 1e4 ({median_total:.2f}) — ensure you are supplying log1p-normalized counts")

    return adata

# ----------------------------
# Main
# ----------------------------
def main(input_file, output_file, cfg, use_gpu=False):
    if input_file.endswith(".h5mu"):
        print("📦 Detected multimodal data (.h5mu) → using RNA modality")
        mdata = mu.read_h5mu(input_file)
        adata = mdata.mod["rna"].copy()
        adata = prepare_adata_for_celltypist(adata)
        adata = run_celltypist_annotation(adata, cfg, use_gpu)
        mdata.mod["rna"] = adata
        mdata.update()
        mdata.write_h5mu(output_file)
    elif input_file.endswith(".h5ad"):
        print("📦 Detected single-modality data (.h5ad)")
        adata = sc.read(input_file)
        adata = prepare_adata_for_celltypist(adata)
        adata = run_celltypist_annotation(adata, cfg, use_gpu)
        adata.write_h5ad(output_file)
    else:
        raise ValueError("Unsupported input format — must be .h5ad or .h5mu")
    print(f"✅ CellTypist annotation completed → {output_file}")


# ----------------------------
# CLI / Snakemake entry point
# ----------------------------
if __name__ == "__main__":
    try:
        # Snakemake execution
        input_file = snakemake.input[0]
        output_file = snakemake.output[0]
        cfg = snakemake.params.config
        use_gpu = snakemake.params.use_gpu
    except NameError:
        # CLI execution
        ap = argparse.ArgumentParser(description="Cell type annotation with CellTypist")
        ap.add_argument("--input", "-i", required=True, help="Input .h5ad or .h5mu file")
        ap.add_argument("--output", "-o", required=True, help="Output annotated file")
        ap.add_argument("--config", "-c", required=True, help="YAML configuration file")
        ap.add_argument("--use_gpu", action="store_true", help="Use GPU (RAPIDS backend)")
        args = ap.parse_args()
        with open(args.config) as f:
            cfg_all = yaml.safe_load(f)
            cfg = cfg_all.get("celltypist", cfg_all)
        input_file, output_file, model_dir, use_gpu = args.input, args.output, args.model_dir, args.use_gpu

    main(input_file, output_file, cfg, use_gpu)
