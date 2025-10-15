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

import scanpy as sc
import muon as mu
import os
import sys
from pathlib import Path
import argparse

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
        print("→ Loaded AnnData (.h5ad).")
        adata.obs["source_file"] = os.path.basename(filepath)
        return adata

    elif ext == ".mtx":
        folder = os.path.dirname(filepath)
        adata = sc.read_10x_mtx(folder, var_names="gene_symbols", cache=True)
        print("→ Loaded 10X matrix (.mtx).")
        adata.obs["source_file"] = os.path.basename(filepath)
        return adata

    elif ext in [".rds", ".h5seurat"]:
        obj = convert_seurat_to_adata(filepath)
        if hasattr(obj, "mod"):
            print("→ Converted Seurat multimodal object to MuData.")
            for mod in obj.mod.values():
                mod.obs["source_file"] = os.path.basename(filepath)
        else:
            print("→ Converted Seurat single-modality object to AnnData.")
            obj.obs["source_file"] = os.path.basename(filepath)
        return obj

    elif ext == ".h5mu":
        mdata = mu.read_h5mu(filepath)
        modalities = list(mdata.mod.keys())
        multimodal = "rna" in modalities and "prot" in modalities
        print(f"→ Loaded MuData (.h5mu) with modalities: {modalities}")
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
    Convert a Seurat .rds or .h5Seurat file into AnnData or MuData
    using MuDataSeurat (R). If multimodal (RNA + ADT), saves to .h5mu;
    if single modality, saves to .h5ad (RNA assay only).

    Returns either an AnnData or MuData object.
    """
    from rpy2.robjects import r
    import tempfile
    import muon as mu

    tmp_dir = tempfile.mkdtemp()
    tmp_h5mu = os.path.join(tmp_dir, "temp.h5mu")
    tmp_h5ad = os.path.join(tmp_dir, "temp.h5ad")

    print("🔁 Converting Seurat object using R/MuDataSeurat...")

    r(f'''
        library(Seurat)
        library(MuDataSeurat)

        # Load Seurat object
        if (grepl("\\.rds$", "{seurat_path}", ignore.case = TRUE)) {{
            cat("Detected .rds file — reading via readRDS()\\n")
            seurat_obj <- readRDS("{seurat_path}")
        }} else {{
            cat("Detected .h5Seurat file — loading via LoadH5Seurat()\\n")
            suppressMessages(library(SeuratDisk))
            seurat_obj <- LoadH5Seurat("{seurat_path}", verbose = FALSE)
        }}

        assays <- names(seurat_obj@assays)
        cat("Detected assays:", assays, "\\n")

        if (all(c("RNA", "ADT") %in% assays)) {{
            cat("Detected multimodal Seurat object (RNA + ADT) — saving as .h5mu\\n")
            # Rename ADT → prot before export for consistency with MuData convention
            names(seurat_obj@assays)[names(seurat_obj@assays) == "ADT"] <- "prot"
            WriteH5MU(seurat_obj, "{tmp_h5mu}")
        }} else if ("RNA" %in% assays) {{
            cat("Single modality (RNA only) — saving as .h5ad\\n")
            WriteH5AD(seurat_obj, "{tmp_h5ad}", assay = "RNA")
        }} else {{
            stop("Unsupported assay combination: only RNA or RNA+ADT supported.")
        }}
    ''')

    # Now load the correct file type in Python
    if os.path.exists(tmp_h5mu):
        print("→ Loaded multimodal MuData (.h5mu)")
        mdata = mu.read_h5mu(tmp_h5mu)
        return mdata
    elif os.path.exists(tmp_h5ad):
        print("→ Loaded single-modality AnnData (.h5ad)")
        adata = sc.read(tmp_h5ad)
        return adata
    else:
        raise FileNotFoundError("Conversion failed — no .h5mu or .h5ad file produced.")


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

    print(f"🧬 Merging {len(adatas)} AnnData objects...")
    merged = adatas[0].concatenate(
        *adatas[1:],
        join="outer",
        batch_key="batch",
        batch_categories=[f"sample_{i+1}" for i in range(len(adatas))]
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
        print(f"🔹 Merging modality: {mod}")
        adatas = [m.mod[mod] for m in mdatas if mod in m.mod]
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
