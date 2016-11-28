FROM ubuntu:16.04
MAINTAINER Peter Humburg <p.humburg@garvan.org.au>

## Install dependencies through distribution packages
RUN apt-get update && apt-get install -y \
    software-properties-common && \
    add-apt-repository ppa:pypy/ppa
RUN apt-get update && apt-get install -y \
    pypy \
    python-pip \
    virtualenv \
    zlib1g-dev \
    zip \
    libncurses5-dev \
    wget \
    git \
    && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

## create virtual environment for pypy
RUN virtualenv -p /usr/bin/pypy /pypy_env && \
    . /pypy_env/bin/activate

## Install required python packages
ADD requirements.txt /tmp/install
RUN pip install -r /tmp/install

## Install aditional packages for debugging and profiling
RUN pip install line_profiler \
                memory_profiler \
                profiling \
                psutil && \
    rm -rf /root/.cache

## Install pyrates package
# using development version mounted from host for now

VOLUME ["/data/fastq", "/data/assembled", "/data/test"]
WORKDIR /code