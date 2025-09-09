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
    libpcre2-dev libpcre3-dev liblzma-dev libbz2-dev libicu-dev zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

# Install Python 3.11 via deadsnakes
RUN add-apt-repository ppa:deadsnakes/ppa \
 && apt-get update && apt-get install -y --no-install-recommends \
    python3.11 python3.11-venv python3.11-dev \
 && curl -sS https://bootstrap.pypa.io/get-pip.py | python3.11 \
 && ln -s /usr/bin/python3.11 /usr/local/bin/python3 \
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
    scikit-learn==1.7.1 \
    plotly==6.3.0 \
    dash==3.2.0 \
    scipy==1.15.3 \
    matplotlib==3.10.6 \
    seaborn==0.13.2 \
    igraph==0.11.9 \
    leidenalg==0.10.2 \
    zarr==3.1.2 \
    loompy==3.0.8 \
    decoupler==2.1.1 \
    gseapy==1.1.9 \
    goatools==1.5.1 \
    pyscenic==0.12.1 \
    celltypist ==1.7.1

# Install Seurat and SeuratDisk (for Seurat -> AnnData conversion)
RUN R -e "install.packages(c('Seurat', 'SeuratObject'), repos='https://cloud.r-project.org')" \
    && R -e "install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('mojaveazure/SeuratDisk')"

# Set working directory
WORKDIR /pipeline

# Default command
ENTRYPOINT ["/bin/bash"]
