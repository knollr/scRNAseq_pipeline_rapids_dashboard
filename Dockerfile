# Use a RAPIDS / CI-Conda image that has Python 3.12 and CUDA 11.4
FROM rapidsai/ci-conda:cuda11.4.3-ubuntu20.04-py3.12  AS base
# (This image exists on Docker Hub) :contentReference[oaicite:0]{index=0}

# Optional: set noninteractive
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
#ENV PATH=/opt/conda/bin:$PATH
ENV CONDA_HTTP_RETRIES=5
ENV CONDA_DOWNLOAD_TIMEOUT=300

# Install extra system dependencies you need
RUN apt-get update && apt-get install -y --no-install-recommends \
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
        libgit2-dev \
        make \
        g++ \
        pkg-config \
        # add any libs you need, e.g. libhdf5, etc.
    && rm -rf /var/lib/apt/lists/*

# Create and activate a conda environment (if not already present)
RUN mkdir -p /opt && \
    wget -qO /tmp/Miniforge3.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash /tmp/Miniforge3.sh -u -b -p /opt/conda && \
    rm /tmp/Miniforge3.sh && \
    /opt/conda/bin/conda clean -afy

RUN /opt/conda/bin/conda create -y -n pipeline python=3.12 \
    && /opt/conda/bin/conda clean -afy

SHELL ["/bin/bash", "-c"]

RUN source /opt/conda/etc/profile.d/conda.sh && conda activate pipeline && \
    conda install -y -c conda-forge mamba && \
    mamba install -y -c rapidsai -c nvidia -c conda-forge -c bioconda \
        rapids=25.06 cudatoolkit=11.4 snakemake && \
    conda clean -afy

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
    R -e "remotes::install_github('mojaveazure/seurat-disk')" && \
    R -e "devtools::install_github('PMBio/MuDataSeurat')"
    

RUN source /opt/conda/etc/profile.d/conda.sh && conda activate pipeline && \
    R -q -e "if ('Seurat' %in% rownames(installed.packages())) { \
                cat('✅ Seurat is installed\n'); \
              } else { \
                cat('❌ Seurat NOT installed\n'); \
                q(status=1); \
              }"


# Now install your additional Python packages
RUN source /opt/conda/etc/profile.d/conda.sh && conda activate pipeline && \
    pip install --no-cache-dir \
        scanpy==1.11.4 anndata==0.12.2 pandas==2.3.2 \
        mudata==0.3.2 muon==0.1.7 hdf5plugin==5.1.0 \
        scikit-learn==1.7.1 plotly==6.3.0 dash==3.2.0 \
        scipy==1.15.3 matplotlib==3.10.6 seaborn==0.13.2 \
        igraph==0.11.9 leidenalg==0.10.2 zarr==3.1.3 \
        loompy==3.0.8 decoupler==2.1.1 gseapy==1.1.9 \
        goatools==1.5.1 pyscenic==0.12.1 celltypist==1.7.1 \
        pyyaml rpy2

# validation of installations
RUN source /opt/conda/etc/profile.d/conda.sh && conda activate pipeline && \
    python -c "import scanpy, anndata, rpy2; print('✅ Python OK')" && \
    R -q -e "library(Seurat); library(SeuratDisk); cat('✅ R OK\n')"

ENV PATH=/opt/conda/envs/pipeline/bin:/opt/conda/bin:$PATH

# Set working directory
WORKDIR /pipeline

# Fully activate conda environment on container start
# Ensure conda.sh is readable and add auto-activation for all shells
RUN chmod -R a+rX /opt/conda/etc/profile.d && \
    echo ". /opt/conda/etc/profile.d/conda.sh && conda activate pipeline" >> /etc/bash.bashrc
#ENTRYPOINT ["/bin/bash"]

# Add convenient shell entry wrapper
RUN echo '#!/bin/bash\n' \
         'source /opt/conda/etc/profile.d/conda.sh\n' \
         'conda activate pipeline\n' \
         'if [ "$1" = "snakemake" ]; then\n' \
         '    shift\n' \
         '    exec snakemake "$@"\n' \
         'else\n' \
         '    exec "$@"\n' \
         'fi' > /usr/local/bin/run_pipeline && chmod +x /usr/local/bin/run_pipeline


RUN mkdir -p /pipeline && chmod -R a+rwx /pipeline
# --- Fix permissions for Singularity users ---
RUN chmod -R a+rX /opt/conda /pipeline /usr/local/bin
RUN ls -ld /opt/conda /opt/conda/bin /opt/conda/envs/pipeline/bin

ENTRYPOINT ["/usr/local/bin/run_pipeline"]
CMD ["/bin/bash"]