# Snakefile

import os
from datetime import datetime
import pathlib

configfile: "config.yaml"


# ----------------------------
# Setup: timestamped results folder
# ----------------------------
run_id = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
results_dir = f"results/run_{run_id}"
os.makedirs(results_dir, exist_ok=True)

# Create or update symlink to latest run
latest_symlink = "results/latest"
if os.path.islink(latest_symlink) or os.path.exists(latest_symlink):
    os.remove(latest_symlink)
# Use absolute paths for cross-directory robustness
os.symlink(os.path.abspath(results_dir), latest_symlink)

# ----------------------------
# Detect file extension for output type
# ----------------------------
ext = os.path.splitext(config["inputs"][0])[1]
ext = (
    ".h5mu"
    if (ext in [".rds", ".h5seurat"] and config.get("multi_modal_seurat", False))
    else (".h5ad" if ext in [".rds", ".h5seurat"] else ext)
)

# ----------------------------
# Rule: Load and merge input files
# ----------------------------
rule load_input:
    """
    Load multiple raw single-cell files and merge into one AnnData or MuData object.
    The output type (.h5ad or .h5mu) is determined automatically by the script.
    """
    input:
        config["inputs"]
    output:
        f"{results_dir}/merged{ext}"
    script:
        "scripts/load_input.py"

# ----------------------------
# Rule: QC and preprocessing
# ----------------------------
rule preprocess:
    """
    Perform QC, normalization, WNN (if multimodal), UMAP, and clustering
    on the merged dataset.
    """
    input:
        f"{results_dir}/merged{ext}"
    output:
        f"{results_dir}/merged_preprocessed{ext}"
    params:
        config=config["preprocessing"],
        use_gpu=config.get("use_gpu", False)
    script:
        "scripts/qc_and_preproc.py"

# ----------------------------
# Rule: automated cell type annotation using CellTypist
# ----------------------------
rule celltypist_annotation:
    """
    Annotate cell types using CellTypist on the preprocessed dataset.
    Works for both single- and multi-modal inputs (.h5ad / .h5mu).
    Downloads CellTypist reference models once and reuses them.
    """
    input:
        f"{results_dir}/merged_preprocessed{ext}"
    output:
        f"{results_dir}/merged_annotated{ext}"
    params:
        config=config["celltypist"],
        use_gpu=config.get("use_gpu", False),
        model_dir="resources/celltypist_models"  
    script:
        "scripts/celltypist_annotation.py"

# ----------------------------
# Rule: all
# ----------------------------
rule all:
    """
    Final target: single preprocessed dataset (h5ad or h5mu)
    """
    input:
        f"{results_dir}/merged_annotated{ext}"
