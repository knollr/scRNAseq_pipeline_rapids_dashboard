# Snakefile

configfile: "config.yaml"

import os
ext = os.path.splitext(config["inputs"][0])[1] 
# attention: currently does not work if input has seurat extension (.rds)! If seurat, output can be both .h5ad or .h5mu depending on modalities!!
# if ext is .rds OR .h5seurat and config["multi_modal_seurat"] is true, then ext = .h5mu else ext = .h5ad, if ext is not .rds, keep ext as is
ext = ".h5mu" if (ext in [".rds", ".h5seurat"] and config.get("multi_modal_seurat", False)) else (".h5ad" if ext in [".rds", ".h5seurat"] else ext)


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
        f"results/merged{ext}"
    script:
        "scripts/load_input.py"
        #shell("python scripts/load_input.py --input {input} --output {output}")

# ----------------------------
# Rule: QC and preprocessing
# ----------------------------
rule preprocess:
    """
    Perform QC, normalization, WNN (if multimodal), UMAP, and clustering
    on the merged dataset.
    """
    input:
        f"results/merged{ext}"
    output:
        f"results/merged_preprocessed{ext}"
    params:
        config=config["preprocessing"]
    script:
        "scripts/qc_and_preproc.py"


# ----------------------------
# Rule: all
# ----------------------------
rule all:
    """
    Final target: single preprocessed dataset (h5ad or h5mu)
    """
    input:
        f"results/merged_preprocessed{ext}"
