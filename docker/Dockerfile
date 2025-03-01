FROM ubuntu:22.04

ENV container=docker
ENV DEBIAN_FRONTEND=noninteractive

ENV SAMTOOLS_VERSION=1.20
ENV PBMM2_VERSION=1.13.1
ENV PBSV_VERSION=2.8.0
ENV PBTK_VERSION=3.4.0
ENV HIPHASE_VERSION=1.4.5
ENV BCFTOOLS_VERSION=1.20
ENV PB_CPG_TOOLS_VERSION=3.0.0
ENV BEDTOOLS_VERSION=2.31.0
ENV HTSLIB_VERSION=1.20
ENV LIBTORCH_VERSION=2.2.0
ENV MOSDEPTH_VERSION=0.3.9

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
    && wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar --no-same-owner -jxvf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2  \
    && cd samtools-${SAMTOOLS_VERSION}/ \
    && ./configure --prefix /opt \
    && make -j && make install

RUN cd /opt/ \
    && wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz \
    && tar --no-same-owner -zxvf bedtools-${BEDTOOLS_VERSION}.tar.gz \
    && rm -rf bedtools-${BEDTOOLS_VERSION}.tar.gz \
    && cd bedtools2 \
    && make

RUN cd /opt/ \
    && wget -O htslib-${HTSLIB_VERSION}.tar.bz2 https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar --no-same-owner -jxvf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && rm htslib-${HTSLIB_VERSION}.tar.bz2 \
    && cd htslib-${HTSLIB_VERSION} \
    && ./configure --prefix /opt \
    && make -j \
    && make install

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

# Other tools
RUN cd /opt \
    && wget https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth \
	&& chmod +x mosdepth

# Install libtorch
RUN cd /opt/ \
    && wget https://download.pytorch.org/libtorch/cu118/libtorch-shared-with-deps-${LIBTORCH_VERSION}%2Bcu118.zip \
    && unzip libtorch-shared-with-deps-${LIBTORCH_VERSION}+cu118.zip \
    && rm libtorch-shared-with-deps-${LIBTORCH_VERSION}+cu118.zip \
    && chmod -R +x libtorch

ENV LIBTORCH_CXX11_ABI=0
ENV LIBTORCH=/opt/libtorch
ENV LD_LIBRARY_PATH=/opt/libtorch/lib
ENV DYLD_LIBRARY_PATH=/opt/libtorch/lib:${LD_LIBRARY_PATH}

# Install Rust and Cargo
ENV RUSTUP_HOME=/opt/rust
ENV CARGO_HOME=/opt/cargo
ENV PATH=${CARGO_HOME}:${LD_LIBRARY_PATH}:${DYLD_LIBRARY_PATH}:${LIBTORCH}${PATH}

RUN curl https://sh.rustup.rs -sSf | sh -s -- -y --no-modify-path \
    && chmod -R a+w ${RUSTUP_HOME} ${CARGO_HOME}
ENV PATH=${CARGO_HOME}/bin/:${PATH}

RUN cd /opt/ \
    && git clone https://github.com/fiberseq/fibertools-rs --branch v0.6.2 \
    && cd fibertools-rs \
    && cargo build --release --all-features \ 
    && mv ./target/release/ft /opt

RUN cd /opt/ \
    && git clone https://github.com/jackh726/bigtools --branch v0.5.5 \
    && cd fibertools-rs \
    && cargo build --release \ 
    && mv ./target/release/* /opt

    ENV PATH=/opt/pbtk/:/opt/samtools-${SAMTOOLS_VERSION}/:/opt/bedtools2/bin:/opt/htslib-${HTSLIB_VERSION}/:/opt/bcftools/:/opt:${PATH}

SHELL ["/bin/bash", "-c"]
CMD ["/bin/bash"]






