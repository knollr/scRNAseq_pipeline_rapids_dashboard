# Use a RAPIDS / CI-Conda image that has Python 3.12 and CUDA 11.4
FROM rapidsai/ci-conda:cuda11.4.3-ubuntu20.04-py3.12  AS base
# (This image exists on Docker Hub) :contentReference[oaicite:0]{index=0}

# Optional: set noninteractive
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=/opt/conda/bin:$PATH
ENV CONDA_HTTP_RETRIES=5
ENV CONDA_DOWNLOAD_TIMEOUT=300

# Install extra system dependencies you need
RUN apt-get update && apt-get install -y --no-install-recommends \
<<<<<<< HEAD
        build-essential \
        git \
        wget \
        curl \
        ca-certificates \
        bzip2 \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libhdf5-dev \
        libfftw3-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        gfortran \
        libblas-dev \
        liblapack-dev \
        build-essential \
        # add any libs you need, e.g. libhdf5, etc.
    && rm -rf /var/lib/apt/lists/*

# Create and activate a conda environment (if not already present)
RUN mkdir -p /opt && \
    wget -qO /tmp/Miniforge3.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash /tmp/Miniforge3.sh -u -b -p /opt/conda && \
    rm /tmp/Miniforge3.sh && \
    /opt/conda/bin/conda clean -afy
=======
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libfftw3-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    gfortran \
    libblas-dev \
    liblapack-dev \
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
>>>>>>> a3bacf1595740b5cd5cf372ab435cb7fdd781e30

RUN /opt/conda/bin/conda create -y -n pipeline python=3.12 \
    && /opt/conda/bin/conda clean -afy

SHELL ["/bin/bash", "-c"]

RUN source /opt/conda/etc/profile.d/conda.sh && conda activate pipeline && \
    conda install -y -c conda-forge mamba && \
    mamba install -y -c rapidsai -c nvidia -c conda-forge \
        rapids=25.06 cudatoolkit=11.4 && \
    conda clean -afy


# Now install your additional Python packages
RUN pip install --no-cache-dir \
<<<<<<< HEAD
        scanpy==1.11.4 anndata==0.12.2 pandas==2.3.2 \
        mudata==0.3.2 hdf5plugin==5.1.0 \
        scikit-learn==1.7.1 plotly==6.3.0 dash==3.2.0 \
        scipy==1.15.3 matplotlib==3.10.6 seaborn==0.13.2 \
        igraph==0.11.9 leidenalg==0.10.2 zarr==3.1.3 \
        loompy==3.0.8 decoupler==2.1.1 gseapy==1.1.9 \
        goatools==1.5.1 pyscenic==0.12.1 celltypist==1.7.1

# Install R + Seurat in the pipeline conda env
RUN source /opt/conda/etc/profile.d/conda.sh && conda activate pipeline && \
    mamba install -y -c conda-forge \
        r-base>=4.3 \
        r-seurat \
        r-remotes \
        r-hdf5r \
        r-devtools \
        r-rcpp \
    && conda clean -afy

# Optional: install seurat-disk from GitHub
RUN source /opt/conda/etc/profile.d/conda.sh && conda activate pipeline && \
    R -e "remotes::install_github('mojaveazure/seurat-disk')"

RUN source /opt/conda/etc/profile.d/conda.sh && conda activate pipeline && \
    R -q -e "if ('Seurat' %in% rownames(installed.packages())) { \
                cat('✅ Seurat is installed\n'); \
              } else { \
                cat('❌ Seurat NOT installed\n'); \
                q(status=1); \
              }"
=======
    scanpy==1.11.4 \
    anndata==0.12.2 \
    pandas==2.3.2 \
    numpy==1.26.4 \
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
    celltypist==1.7.1

# Install Seurat and packages
RUN R -e "install.packages('Seurat', repos='https://cloud.r-project.org')"
RUN R -q -e "stopifnot('Seurat' %in% rownames(installed.packages()))"
RUN R -e "install.packages('SeuratObject', repos='https://cloud.r-project.org')"
RUN R -q -e "stopifnot('SeuratObject' %in% rownames(installed.packages()))"
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('mojaveazure/seurat-disk')"
RUN R -q -e "stopifnot('remotes' %in% rownames(installed.packages()))"
RUN R -q -e "stopifnot('SeuratDisk' %in% rownames(installed.packages()))"
>>>>>>> a3bacf1595740b5cd5cf372ab435cb7fdd781e30

# Set working directory
WORKDIR /pipeline

# Fully activate conda environment on container start
ENTRYPOINT ["/bin/bash", "-c", "source /opt/conda/etc/profile.d/conda.sh && conda activate pipeline && exec bash"]