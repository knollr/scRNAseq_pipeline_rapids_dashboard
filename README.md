# Pyhton scRNA-seq Pipeline + Plotly Dashboard

An end-to-end single-cell RNA-seq workflow using:
- **Snakemake** for reproducibility
- **Scanpy + RAPIDS (GPU acceleartion)** for analysis, including:
  - normalization
  - dimensionality reduction
  - ...
- **Plotly Dash** for interactive exploration

As input, either a demultiplexed seurat or adata object can be used.
