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
    software-properties-common \
    r-base \
 && add-apt-repository ppa:deadsnakes/ppa \
 && apt-get update && apt-get install -y --no-install-recommends \
    python3.12 python3.12-venv python3.12-dev \
 && curl -sS https://bootstrap.pypa.io/get-pip.py | python3.12 \
 && ln -s /usr/bin/python3.12 /usr/local/bin/python3 \
 && python3 -m pip install --upgrade pip setuptools wheel \
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
    rmm-cu11==25.6.0 \
    rapids-singlecell==0.13.2 \
    rpy2==3.6.3 \
    scikit-learn==1.3.2 \
    plotly==5.15.0 \
    dash==2.19.0 \
    numpy==1.26.4 \
    scipy==1.13.1 \
    matplotlib==3.9.2 \
    seaborn==0.13.2 \
    scikit-learn==1.5.2 \
    igraph==0.11.8 \
    leidenalg==0.10.2 \
    plotly==5.24.1 \
    zarr==2.18.2 \
    loompy==3.0.7 \
    decoupler==1.6.0 \
    gseapy==1.1.4 \
    goatools==1.4.12 \
    pyscenic==0.12.1

# Install Seurat and SeuratDisk (for Seurat -> AnnData conversion)
RUN R -e "install.packages(c('Seurat', 'SeuratObject'), repos='https://cloud.r-project.org')" \
    && R -e "install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('mojaveazure/SeuratDisk')"

# Set working directory
WORKDIR /pipeline

# Default command
ENTRYPOINT ["/bin/bash"]
