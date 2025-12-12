# 🧬 scRNAseq Pipeline + Dashboard 🚀

A **fast, GPU-accelerated single-cell RNA-seq & multi-omics analysis pipeline** with an interactive dashboard for easy exploration. Built with **Snakemake**, **Rapids**, **Scanpy**, and **MuData**. Perfect for researchers who want speed without sacrificing flexibility!  

### ⚡ Features
- 🖥️ **GPU-optimized** processing with automatic fallback to CPU  
- 🔄 **Seamless data conversion**: Seurat (`.rds` / `.h5Seurat`) → AnnData / MuData
- 🔗 **Merging of multiple datasets** for integrated analysis  
- 🔬 **Multi-omics support**: RNA and protein (CITE/ADT)
- 🧪 **Automated cell type annotation**: Cell type annotation using CellTypist 
- 📊 **Interactive dashboard** for visual exploration and QC 
- 🧩 Modular **Snakemake workflow** for reproducibility  

### 🚀 Quick Start
1. 🐳 Use the provided **Docker container** for zero-hassle setup ([knollr/scrnaseq_pipeline:latest](https://hub.docker.com/r/knollr/scrnaseq_pipeline/tags))
2. 📁 Get your dataset(s), different formats possible (`.h5ad`, `.h5mu`, `.rds`, `.h5seurat`)
3. 📝 Set your configuration file (`config.yaml`) with dataset paths and parameters  
4. ⚡ Run the snakemake pipeline and explore your results in the dashboard  

### 🔧 Tech Stack
- Python and R integration  
- Scanpy, MuData (incl dataset merging)  
- Rapids Single Cell / GPU-accelerated clustering
- Snakemake workflow manager
