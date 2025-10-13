# Snakefile

configfile: "config.yaml"

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
        # Temporary placeholder; will be replaced by actual path in 'preprocess' via 'expand'
        "results/merged/merged"
    run:
        # Call Python script
        shell("python scripts/load_input.py --input {input} --output {output}")

        # Detect whether the script saved .h5ad or .h5mu
        import os
        h5ad_path = f"{output}.h5ad"
        h5mu_path = f"{output}.h5mu"
        if os.path.exists(h5ad_path):
            # Rename to final output path
            shell(f"mv {h5ad_path} {output}.h5ad")
            touch(output + ".type")  # create a marker file
        elif os.path.exists(h5mu_path):
            shell(f"mv {h5mu_path} {output}.h5mu")
            touch(output + ".type")
        else:
            raise ValueError("load_input.py did not produce a .h5ad or .h5mu file")

# ----------------------------
# Rule: QC and preprocessing
# ----------------------------
rule preprocess:
    """
    Perform QC, normalization, WNN (if multimodal), UMAP, and clustering
    on the merged dataset.
    """
    input:
        merged=lambda wildcards: 
            # find merged file dynamically
            next(
                f for f in ["results/merged/merged.h5ad", "results/merged/merged.h5mu"]
                if os.path.exists(f)
            )
    output:
        preprocessed=lambda wildcards:
            input.merged.replace("merged", "merged_preprocessed")
    params:
        config=config["preprocessing"]
    shell:
        """
        python scripts/qc_and_preproc.py \
            --input {input.merged} \
            --output {output.preprocessed} \
            --config {params.config}
        """

# ----------------------------
# Rule: all
# ----------------------------
rule all:
    """
    Final target: single preprocessed dataset (h5ad or h5mu)
    """
    input:
        lambda wildcards: next(
            f for f in ["results/preprocessed/merged_preprocessed.h5ad",
                        "results/preprocessed/merged_preprocessed.h5mu"]
            if os.path.exists(f)
        )
