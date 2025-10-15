# Snakefile

configfile: "config.yaml"

import os
ext = os.path.splitext(config["inputs"][0])[1] 

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
