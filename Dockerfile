# Base image with CUDA 11.4 for GPU
FROM nvidia/cuda:11.4.3-cudnn8-runtime-ubuntu20.04

# Set environment
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    git wget curl bzip2 ca-certificates \
    build-essential cmake pkg-config \
    libhdf5-serial-dev hdf5-tools \
    python3 python3-pip python3-dev python3-venv \
    r-base \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip
RUN python3 -m pip install --upgrade pip setuptools wheel

# Install Python packages with fixed versions (recent stable)
RUN pip install --no-cache-dir \
    scanpy==1.11.4 \
    anndata==0.12.2 \
    pandas==2.3.2 \
    numpy==2.3.2 \
    mudata==0.3.2 \
    hdf5plugin==5.1.0 \
    rmm==23.12 \
    rapids-singlecell==0.13.2 \
    rpy2==3.7.0 \
    scikit-learn==1.3.2 \
    plotly==5.15.0 \
    dash==2.19.0

# Install SeuratDisk (to convert Seurat .rds to AnnData)
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('mojaveazure/SeuratDisk')"

# Set working directory
WORKDIR /pipeline

# Default command
ENTRYPOINT ["/bin/bash"]
