FROM --platform=linux/amd64 ubuntu:20.04

# Set environment variables for non-interactive installs, path, and python
ENV DEBIAN_FRONTEND=noninteractive \
    JAVA_HOME=/opt/java/jdk18 \
    MINICONDA_VERSION=py312_24.1.2-0 \
    PATH=/opt/conda/envs/SpaceTracer/bin:/opt/conda/bin:/opt/java/jdk18/bin:$PATH

# Install base dependencies, JDK, and Miniconda in a single layer to reduce image size
RUN apt-get update && \
    apt-get install -y --no-install-recommends wget bzip2 build-essential openssh-server curl ca-certificates zlib1g-dev locales jq && \
    locale-gen en_US.UTF-8 && \
    # Install JDK
    mkdir -p /opt/java && \
    cd /opt/java && \
    curl -L https://download.java.net/java/GA/jdk18.0.1.1/65ae32619e2f40f3a9af3af1851d6e19/2/GPL/openjdk-18.0.1.1_linux-x64_bin.tar.gz -o jdk18.tar.gz && \
    tar -xzf jdk18.tar.gz && \
    rm jdk18.tar.gz && \
    mv jdk-18.0.1.1 jdk18 && \
    # Install Miniconda
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    # Initialize conda & setup channels/solvers
    /opt/conda/bin/conda clean --all --yes && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda config --add channels anaconda && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels r && \
    # Install Mamba to speed up future solves
    /opt/conda/bin/conda install -y -n base -c conda-forge mamba && \
    # Clean up APT cache
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Create environments and install packages with mamba, then aggressively clean cache
RUN mamba create -y -n snakemake3 -c conda-forge -c bioconda -c defaults snakemake=3.13.3 && \
    mamba create -y -n SpaceTracer python=3.9 && \
    mamba install -y -n SpaceTracer numpy=1.24.3 pandas=1.5.3 scikit-learn matplotlib=3.5.0 scipy=1.11.3 seaborn=0.13.0 && \
    mamba install -y -n SpaceTracer samtools bedtools pybedtools vcftools geopandas=0.14.0 pip=21.2.4 && \
    mamba install -y -n SpaceTracer r-base r-pracma r-dplyr r-deconstructSigs r-extraDistr r-tidyr && \
    mamba clean --all --yes

# Install pip dependencies within the SpaceTracer environment
# Using --no-cache-dir reduces final image size
RUN /bin/bash -c "source activate SpaceTracer && \
    pip install --no-cache-dir uv && \
    pip install --no-cache-dir --no-build-isolation umi-tools==1.1.6 && \
    pip install --no-cache-dir sinto && \
    uv pip install --no-cache colorcet==3.0.1 \
    GraphST hyperopt imblearn imbalanced_learn SpaGCN==1.2.7 \
    vcfpy==0.13.6 torch esda==2.5.0 wordcloud joblib==1.3.2 \
    opencv_python pyfaidx==0.7.2.2 pysal==23.07 pysam==0.22.0 \
    scanpy==1.10.2 tqdm==4.66.1 shap"

# Setup workflow directory and SSH Daemon
RUN mkdir /mnt/workflow && \
    chown root:root /mnt/workflow && \
    chmod 755 /mnt/workflow && \
    mkdir -p /var/run/sshd && \
    echo "source activate SpaceTracer" >> ~/.bashrc

WORKDIR /mnt/workflow
SHELL ["/bin/bash", "-c"]
CMD ["/usr/sbin/sshd", "-D"]
