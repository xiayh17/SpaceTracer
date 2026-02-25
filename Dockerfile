FROM ubuntu:20.04

# Set DEBIAN_FRONTEND to noninteractive to avoid tzdata prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install wget, OpenSSH, and other dependencies
RUN apt-get update && apt-get install -y wget bzip2 build-essential openssh-server curl zlib1g-dev locales jq && locale-gen en_US.UTF-8 && rm -rf /var/lib/apt/lists/*

# Install JDK 18.0.1.1
RUN mkdir -p /opt/java && \
    cd /opt/java && \
    curl -L https://download.java.net/java/GA/jdk18.0.1.1/65ae32619e2f40f3a9af3af1851d6e19/2/GPL/openjdk-18.0.1.1_linux-x64_bin.tar.gz -o jdk18.tar.gz && \
    tar -xzf jdk18.tar.gz && \
    rm jdk18.tar.gz && \
    mv jdk-18.0.1.1 jdk18

# Set JAVA_HOME and add to PATH
ENV JAVA_HOME=/opt/java/jdk18
ENV PATH=$JAVA_HOME/bin:$PATH

# Install miniconda
ENV MINICONDA_VERSION=py312_24.1.2-0
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -t -p -y

# Add conda to PATH
ENV PATH=/opt/conda/bin:$PATH

# conda init
RUN conda init bash

# Config conda with libmamba solver for faster installation
RUN conda update -n base -c defaults conda && \
    conda install -y -n base -c conda-forge conda-libmamba-solver && \
    conda config --set solver libmamba

# Create snakemake environment and install snakemake3.13.3
RUN conda create -y -n snakemake3 -c bioconda snakemake=3.13.3

# Create SpaceTracer environment and configure channels
RUN conda config --add channels conda-forge && \
    conda config --add channels anaconda && \
    conda config --add channels bioconda && \
    conda config --add channels r && \
    conda create -y -n SpaceTracer python=3.9

# Install conda packages with libmamba solver for faster installation
RUN conda install -y -n SpaceTracer \
    numpy=1.24.3 \
    pandas=1.5.3 \
    scikit-learn \
    matplotlib=3.5.0 \
    scipy=1.11.3 \
    seaborn=0.13.0 \
    samtools \
    bedtools \
    pybedtools \
    vcftools \
    geopandas=0.14.0 \
    r-base \
    pip=21.2.4 \
    r-pracma \
    r-dplyr \
    r-deconstructSigs \
    r-extraDistr \
    r-tidyr && \
    conda clean --all --yes

# Install uv for faster pip package installation
RUN pip install uv

# Install pip packages with uv for faster installation
RUN /bin/bash -c "source activate SpaceTracer && \
    pip install --no-build-isolation umi-tools==1.1.6 && \
    pip install sinto && \
    uv pip install colorcet==3.0.1 \
    GraphST \
    hyperopt \
    imblearn \
    imbalanced_learn \
    SpaGCN==1.2.7 \
    vcfpy==0.13.6 \
    torch \
    esda==2.5.0 \
    wordcloud \
    joblib==1.3.2 \
    opencv_python \
    pyfaidx==0.7.2.2 \
    pysal==23.07 \
    pysam==0.22.0 \
    scanpy==1.10.2 \
    tqdm==4.66.1 \
    shap"

# Create a directory for mounting the workflow
RUN mkdir /mnt/workflow && \
    chown root:root /mnt/workflow && \
    chmod 755 /mnt/workflow

WORKDIR /mnt/workflow

# Set the conda environment to be activated by default
SHELL ["/bin/bash", "-c"]
RUN echo "conda activate SpaceTracer" >> ~/.bashrc
ENV PATH=/opt/conda/envs/SpaceTracer/bin:$PATH

CMD ["/usr/sbin/sshd", "-D"]
