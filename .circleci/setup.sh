#!/bin/bash

go version

# Install Singularity
export VERSION=3.5.3 && # adjust this as necessary \
    wget https://dabdceba-6d04-11e5-ba46-22000b92c6ec.e.globus.org/containers/public/singularity-${VERSION}.tar.gz && \
    tar -xzf singularity-${VERSION}.tar.gz && \
    cd singularity/
./mconfig && \
    make -C ./builddir && \
    sudo make -C ./builddir install

singularity --version
