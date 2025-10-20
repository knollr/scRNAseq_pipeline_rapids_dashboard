#!/usr/bin/env python3
"""
Load one or multiple single-cell objects (AnnData, MuData, or Seurat)
and merge them into a single AnnData or MuData object without performing QC.

Supports:
- .h5ad       : AnnData
- .mtx        : 10X Genomics format
- .rds/.h5Seurat : Seurat object (converted via SeuratDisk)
- .h5mu       : MuData (multi-modal)
"""
#import scripts.silence_warnings

import os
import sys
import scanpy as sc
import muon as mu
from pathlib import Path
import argparse
import anndata as ad

# ----------------------------
# Functions
# ----------------------------

def load_single_file(filepath: str):
    """
    Load a single file into an AnnData or MuData object.
    Supports .h5ad, .mtx, .rds, .h5Seurat, and .h5mu formats.
    """
    _, ext = os.path.splitext(filepath)
    ext = ext.lower()
    print(f"📂 Loading file: {filepath}")

    if ext == ".h5ad":
        adata = sc.read(filepath)
        print("📂 Loaded AnnData (.h5ad).")
        adata.obs["source_file"] = os.path.basename(filepath)
        return adata

    elif ext == ".mtx":
        folder = os.path.dirname(filepath)
        adata = sc.read_10x_mtx(folder, var_names="gene_symbols", cache=True)
        print("📂 Loaded 10X matrix (.mtx).")
        adata.obs["source_file"] = os.path.basename(filepath)
        return adata

    elif ext in [".rds", ".h5seurat"]:
        obj = convert_seurat_to_adata(filepath)
        if hasattr(obj, "mod"):
            print("🔁 Converted Seurat multimodal object to MuData.")
            for mod in obj.mod.values():
                mod.obs["source_file"] = os.path.basename(filepath)
        else:
            print("🔁 Converted Seurat single-modality object to AnnData.")
            obj.obs["source_file"] = os.path.basename(filepath)
        return obj

    elif ext == ".h5mu":
        mdata = mu.read_h5mu(filepath)
        modalities = list(mdata.mod.keys())
        multimodal = "rna" in modalities and "prot" in modalities
        print(f"📂 Loaded MuData (.h5mu) with modalities: {modalities}")
        if multimodal:
            print("✅ Detected multimodal data: RNA + Protein.")
        else:
            print("⚠️ Single modality MuData detected.")
        # attach filename for each modality
        for mod in mdata.mod.values():
            mod.obs["source_file"] = os.path.basename(filepath)
        return mdata

    else:
        raise ValueError(f"Unsupported input file format: {ext}")


def convert_seurat_to_adata(seurat_path: str):
    """
    Convert a Seurat .rds or .h5Seurat file into AnnData or MuData.
    Handles both Assay (v4) and Assay5 (v5) objects.
    For Seurat v5 multimodal (RNA + ADT), exports via zellkonverter and merges in Python.
    """

    from rpy2.robjects import r
    import tempfile
    import muon as mu
    import scanpy as sc
    from rpy2.rinterface_lib import callbacks
    import sys
    import anndata as ad

    # Capture R console output
    def console_to_stdout(x):
        sys.stdout.write(x)
        sys.stdout.flush()
    callbacks.consolewrite_print = console_to_stdout
    callbacks.consolewrite_warnerror = console_to_stdout

    tmp_dir = tempfile.mkdtemp()
    tmp_h5ad = os.path.join(tmp_dir, "temp.h5ad")
    tmp_h5mu = os.path.join(tmp_dir, "temp.h5mu")
    tmp_h5ad_rna = os.path.join(tmp_dir, "rna_tmp.h5ad")
    tmp_h5ad_prot = os.path.join(tmp_dir, "prot_tmp.h5ad")

    print(f"📁 Temporary directory: {tmp_dir}")
    print("🔁 Converting Seurat object using R (MuDataSeurat / zellkonverter)...")

    # R code block
    r(f"""
        suppressMessages(library(Seurat))
        suppressMessages(library(MuDataSeurat))
        suppressMessages(library(zellkonverter))
        suppressMessages(library(SingleCellExperiment))
        suppressMessages(library(SeuratDisk))
        suppressMessages(library(reticulate))
                
        #use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
        
        ## careful, add at end of dockerfile
        # ENV BASILISK_USE_SYSTEM_DIR=1
        # ENV BASILISK_EXTERNAL_CONDA=1
        # ENV RETICULATE_PYTHON=/opt/conda/envs/pipeline/bin/python

        if (grepl("\\\\.rds$", "{seurat_path}", ignore.case = TRUE)) {{
            message("🔍 Detected .rds file — reading via readRDS()")
            seurat_obj <- readRDS("{seurat_path}")
        }} else {{
            message("🔍 Detected .h5Seurat file — loading via LoadH5Seurat()")
            seurat_obj <- LoadH5Seurat("{seurat_path}", verbose = FALSE)
        }}

        assays <- names(seurat_obj@assays)
        message(paste("🔍 Detected assays:", paste(assays, collapse = ", ")))
        assay_classes <- sapply(seurat_obj@assays, class)
        message(paste("🧩 Assay classes:", paste(unique(assay_classes), collapse = ", ")))

        write_h5ad_safe <- function(obj, filename, assay = "RNA") {{
            message(paste("💾 Writing", filename))
            sce <- as.SingleCellExperiment(obj)
            zellkonverter::writeH5AD(sce, filename)
        }}
        
        # write_h5ad_safe <- function(obj, filename, assay = "RNA") {{
        #     message(paste("💾 Writing", filename, "via reticulate + anndata"))
        #     sce <- as.SingleCellExperiment(obj)
        #     ad <- import("anndata")
        #     np <- import("numpy")
        #     pd <- import("pandas")

        #     # Extract matrix and metadata
        #     X <- assay(sce)
        #     obs <- as.data.frame(colData(sce))
        #     var <- as.data.frame(rowData(sce))

        #     # Convert to numpy arrays / DataFrames
        #     adata <- ad$AnnData(
        #         X = np$array(X),
        #         obs = pd$DataFrame(obs),
        #         var = pd$DataFrame(var)
        #     )
        #     adata$write_h5ad(filename)
        # }}

        if (all(c("RNA", "ADT") %in% assays)) {{
            message("🔬🔍 Multimodal (RNA + ADT) detected")
            names(seurat_obj@assays)[names(seurat_obj@assays) == "ADT"] <- "prot"
            if (any(assay_classes == "Assay5")) {{
                message("⚙️🔍 Assay5 detected — exporting each modality via zellkonverter")
                write_h5ad_safe(seurat_obj[["RNA"]], "{tmp_h5ad_rna}")
                write_h5ad_safe(seurat_obj[["prot"]], "{tmp_h5ad_prot}")
            }} else {{
                message("✅ Using MuDataSeurat::WriteH5MU() for Seurat v4-style object")
                MuDataSeurat::WriteH5MU(seurat_obj, "{tmp_h5mu}")
            }}
        }} else if ("RNA" %in% assays) {{
            message("🧬🔍 single modality (RNA only) detected")
            if (any(assay_classes == "Assay5")) {{
                message("⚙️🔍 Assay5 detected — exporting via zellkonverter")
                write_h5ad_safe(seurat_obj, "{tmp_h5ad}")
            }} else {{
                message("✅ Using MuDataSeurat::WriteH5AD() for Seurat v4-style object")
                MuDataSeurat::WriteH5AD(seurat_obj, "{tmp_h5ad}", assay = "RNA")
            }}
        }} else {{
            stop("❌ Unsupported assay combination: only RNA or RNA+ADT supported.")
        }}
    """)

    # -----------------------------
    # Back to Python: load and merge if multimodal
    # -----------------------------
    if os.path.exists(tmp_h5mu):
        print("📂 Loaded multimodal MuData (.h5mu)")
        mdata = mu.read_h5mu(tmp_h5mu)

        for key in list(mdata.mod.keys()):  # use list() to avoid runtime dict size change
            new_key = key.lower()
            if new_key != key:
                mdata.mod[new_key] = mdata.mod.pop(key)
        print("✅ Created MuData object.")
        return mdata

    elif os.path.exists(tmp_h5ad):
        print("📂 Loaded single-modality AnnData (.h5ad)")
        adata = sc.read(tmp_h5ad)
        print("✅ Created AnnData object.")
        return adata

    elif os.path.exists(tmp_h5ad_rna) and os.path.exists(tmp_h5ad_prot):
        print("🔗 Merging RNA + Protein modalities into MuData (.h5mu)...")
        rna = sc.read(tmp_h5ad_rna)
        prot = sc.read(tmp_h5ad_prot)
        mdata = mu.MuData({"rna": rna, "prot": prot})
        print("✅ Created merged MuData object in memory.")
        return mdata

    else:
        raise FileNotFoundError("❌Conversion failed — no output file produced.")


def merge_adatas(adatas):
    """
    Merge multiple AnnData objects. Performs outer join on variables.
    If MuData objects are detected, merges modalities separately.
    """
    if len(adatas) == 1:
        return adatas[0]

    # Check if any MuData objects are in the list
    if any(hasattr(a, "mod") for a in adatas):
        print("🧩 Detected MuData objects — attempting multimodal merge.")
        merged = merge_mudatas(adatas)
        return merged

    print(f"🔗 Merging {len(adatas)} AnnData objects...")
    merged = ad.concat(
        adatas,
        join="outer",
        label="batch",
        keys=[f"sample_{i+1}" for i in range(len(adatas))]
    )
    print(f"✅ Merged AnnData: {merged.n_obs} cells × {merged.n_vars} genes")
    return merged


def merge_mudatas(mdatas):
    """
    Merge multiple MuData objects by merging shared modalities (e.g. RNA, prot).
    """
    # Collect modality names
    all_modalities = set()
    for m in mdatas:
        all_modalities.update(m.mod.keys())

    merged_mods = {}
    for mod in all_modalities:
        print(f"🔗 Merging modality: {mod}")
        adatas = [m.mod[mod] for m in mdatas if mod in m.mod and not hasattr(m.mod[mod], "mod")]
        merged_mods[mod] = merge_adatas(adatas)

    merged_mdata = mu.MuData(merged_mods)
    print(f"✅ Merged MuData with modalities: {list(merged_mdata.mod.keys())}")
    return merged_mdata

# ----------------------------
# Main
# ----------------------------

def main(input_files, output_file):
    """
    Main function: load and merge input files into one AnnData or MuData object.
    """
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)

    loaded_objects = [load_single_file(f) for f in input_files]
    merged = merge_adatas(loaded_objects)

    # Choose correct write method
    # Determine output path with proper extension
    if hasattr(merged, "write_h5mu"):
        #output_file = output_file + ".h5mu"
        merged.write_h5mu(output_file)
    else:
        #output_file = output_file + ".h5ad"
        merged.write(output_file)

    print(f"💾 Saved merged object to: {output_file}")

# ----------------------------
# CLI / Snakemake entry point
# ----------------------------

if __name__ == "__main__":
    # Handle both Snakemake and CLI usage
    try:
        input_files = list(snakemake.input)
        output_file = snakemake.output[0]
    except NameError:
        p = argparse.ArgumentParser(description="Load and merge single-cell files.")
        p.add_argument("--input", "-i", nargs="+", required=True, help="Input file(s)")
        p.add_argument("--output", "-o", required=True, help="Output .h5ad or .h5mu file")
        args = p.parse_args()
        input_files = args.input
        output_file = args.output

    main(input_files, output_file)
