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

ENV PBMM2_VERSION=1.13.1
ENV PBSV_VERSION=2.9.0
ENV PBTK_VERSION=3.4.0
ENV HIPHASE_VERSION=1.4.5
ENV BCFTOOLS_VERSION=1.20
ENV PB_CPG_TOOLS_VERSION=3.0.0


# Install dependencies
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
    nano \
    python3 python3-pip python-is-python3 \
    jq \
    gawk \
    parallel \
    pkg-config \
    openjdk-21-jre \
    && rm -rf /var/lib/apt/lists/*


RUN cd /opt/ \
    && wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && tar --no-same-owner -jxvf bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && rm bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && cd bcftools-${BCFTOOLS_VERSION} \
    && ./configure --prefix /opt \
    && make -j \
    && make install
    
# Install pacbio tools
RUN wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v${PB_CPG_TOOLS_VERSION}/pb-CpG-tools-v${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
	tar --no-same-owner -zxvf pb-CpG-tools-v${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt && \
	rm pb-CpG-tools-v${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
    mv /opt/pb-CpG-tools-v${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores /opt

RUN wget https://github.com/PacificBiosciences/pbsv/releases/download/v${PBSV_VERSION}/pbsv --directory-prefix /opt \
    && chmod +x /opt/pbsv

RUN mkdir /opt/pbtk/ && cd /opt/pbtk/ \
    && wget https://github.com/PacificBiosciences/pbtk/releases/download/v${PBTK_VERSION}/pbtk.tar.gz \
    && tar --no-same-owner -zvxf pbtk.tar.gz \
    && rm pbtk.tar.gz \
    && chmod +x *

RUN wget https://github.com/PacificBiosciences/pbmm2/releases/download/v${PBMM2_VERSION}/pbmm2 --directory-prefix /opt \
    && chmod +x /opt/pbmm2

RUN wget https://github.com/PacificBiosciences/HiPhase/releases/download/v${HIPHASE_VERSION}/hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
	tar --no-same-owner -zxvf hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt && \
	rm hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
    mv /opt/hiphase-v${HIPHASE_VERSION}-x86_64-unknown-linux-gnu/hiphase /opt


ENV PATH=/opt/bcftools/:/opt:/venv/bin:${PATH}



SHELL ["/bin/bash", "-c"]
CMD ["/bin/bash"]






