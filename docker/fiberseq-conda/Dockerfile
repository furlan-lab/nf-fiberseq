FROM condaforge/miniforge3:latest AS build 
ENV container=docker

ENV DEBIAN_FRONTEND=noninteractive

RUN conda config --set channel_alias https://conda-forge.fredhutch.org \
    && conda config --remove-key channels \
    && conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda config --remove channels defaults \
    && conda config --set channel_priority flexible

SHELL ["/bin/bash", "-c"]

RUN conda create --name fiberseq \
    && conda init bash \
    && source activate base \
    && conda activate fiberseq

RUN conda install -n fiberseq samtools==1.20 
RUN conda install -n fiberseq htslib==1.20 
RUN conda install -n fiberseq bedtools==2.31.0 
RUN conda install -n fiberseq mosdepth
RUN conda install -n fiberseq bioawk ripgrep csvtk 
RUN conda install -n fiberseq hck>=0.9.2
RUN conda install -n fiberseq fibertools-rs==0.6
RUN conda install -n fiberseq bigtools==0.5.4
RUN conda install -c conda-forge conda-pack

RUN conda-pack -n fiberseq -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

RUN /venv/bin/conda-unpack

FROM debian:latest AS runtime

COPY --from=build /venv /venv

ENV BCFTOOLS_VERSION=1.20
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    ca-certificates \
    zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-openssl-dev libssl-dev libncurses5-dev libncursesw5-dev \
    autoconf automake make cmake \
    gcc g++ \
    bzip2 zip unzip tar xz-utils \
    parallel \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

RUN cd /opt/ \
    && wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && tar --no-same-owner -jxvf bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && rm bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && cd bcftools-${BCFTOOLS_VERSION} \
    && ./configure --prefix /opt \
    && make -j \
    && make install

ENV PATH=/opt/bin:/venv/bin:$PATH
CMD ["/bin/bash"]


