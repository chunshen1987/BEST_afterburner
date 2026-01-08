FROM ubuntu:22.04

MAINTAINER Chun Shen <chunshen@wayne.edu>

RUN export DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true && \
    apt-get update && apt-get install software-properties-common -y \
    && add-apt-repository ppa:ubuntu-toolchain-r/test

RUN export DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true && \
    apt-get update && apt-get install -y --no-install-recommends \
        cmake \
        curl \
        g++-13 \
        gcc-13 \
        gfortran-13 \
        git \
        hdf5-tools \
        libfftw3-dev \
        libgsl-dev \
        libhdf5-dev \
        libeigen3-dev \
        libz-dev \
        make \
        pkg-config \
        python3 \
        python3-pip \
        time \
        vim \
        wget \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install pre-reqs python packages
RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install numpy h5py scipy

ENV HOME /opt/

# download and compile the iEBE-MUSIC framework
RUN alias g++=g++-13 && alias gcc=gcc-13 \
&& cd ${HOME} \
&& git clone --depth=1 https://github.com/chunshen1987/BEST_afterburner -b SMASH-3.2.2 \
&& cd BEST_afterburner/external_codes \
&& bash get_iSS.sh \
&& bash get_microcanonical_sampler.sh \
&& bash get_msu_sampler.sh \
&& bash get_smash.sh 

ENV PYTHONIOENCODING utf-8

WORKDIR ${HOME}/BEST_afterburner

CMD /bin/bash
