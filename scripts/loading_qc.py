import scanpy as sc
from pathlib import Path
import sys
import os

# Try Snakemake placeholders
try:
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    qc_params = snakemake.config["qc"]
except NameError:
    # For testing standalone
    input_file = sys.argv[1] if len(sys.argv) > 1 else "data/raw/10x_matrix.mtx"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "results/qc/adata_qc.h5ad"
    qc_params = {"min_genes": 200, "max_genes": 6000, "max_mt": 0.1}

# Detect file type
_, ext = os.path.splitext(input_file)

if ext == ".h5ad":
    adata = sc.read(input_file)
    print("Loaded AnnData object.")
elif ext == ".mtx":
    # Expect accompanying barcodes.tsv and features.tsv in the same folder
    folder = os.path.dirname(input_file)
    adata = sc.read_10x_mtx(folder, var_names='gene_symbols', cache=True)
    print("Loaded 10X matrix.")
elif ext == ".rds":
    import rpy2.robjects as robjects

    print("Converting Seurat .rds object to AnnData...")

    # Define temp filenames
    tmp_h5seurat = output_file.replace(".h5ad", ".h5Seurat")
    tmp_h5ad = output_file

    # Run Seurat -> h5Seurat -> h5ad conversion in R
    robjects.r(f'''
    #library(Seurat)
    library(SeuratDisk)
    seurat_obj <- readRDS("{input_file}")
    SaveH5Seurat(seurat_obj, filename="{tmp_h5seurat}", overwrite=TRUE)
    Convert("{tmp_h5seurat}", dest = "h5ad", overwrite=TRUE)
    ''')

    # Load into AnnData
    adata = sc.read(tmp_h5ad)
    print("Seurat .rds successfully converted to AnnData.")
else:
    raise ValueError(f"Unsupported input file type: {ext}")

# Flag mitochondrial genes
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

# Calculate QC metrics including mt
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# Filter cells
adata = adata[adata.obs.n_genes_by_counts >= qc_params["min_genes"], :]
adata = adata[adata.obs.n_genes_by_counts <= qc_params["max_genes"], :]
adata = adata[adata.obs.pct_counts_mt <= qc_params["max_mt"], :]

# Save processed AnnData
Path(output_file).parent.mkdir(parents=True, exist_ok=True)
adata.write(output_file)
print(f"QC complete: {adata.n_obs} cells, {adata.n_vars} genes saved to {output_file}")