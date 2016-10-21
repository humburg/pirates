FROM ubuntu:16.04
MAINTAINER Peter Humburg <p.humburg@garvan.org.au>

## Install dependencies through distribution packages
RUN apt-get update && apt-get install -y \
    python \
    python-nose \
    python-pip \
    zlib1g-dev \
    zip \
    libncurses5-dev \
    wget \
    git \
    && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

## Install latest version of python packages
RUN pip install \
    argparse \
    coverage \
    versioneer

## Install pyrates package
# using development version mounted from host for now

VOLUME ['/data/fastq', '/data/assembled', '/data/test']
WORKDIR /code