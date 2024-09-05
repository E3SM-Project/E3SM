#!/bin/bash

go version

export VERSION=3.5.3

function download
{
    wget -t 3 -O singularity-${VERSION}.tar.gz http://portal.nersc.gov/project/e3sm/lukasz/singularity-${VERSION}.tar.gz || \
        wget -t 3 -O singularity-${VERSION}.tar.gz https://dabdceba-6d04-11e5-ba46-22000b92c6ec.e.globus.org/containers/public/singularity-${VERSION}.tar.gz
    return $?
}

# Install Singularity deps:

sudo apt-get update && sudo apt-get -y install uuid-dev

# Install Singularity
download && \
    tar -xzf singularity-${VERSION}.tar.gz && \
    cd singularity/
./mconfig && \
    make -C ./builddir && \
    sudo make -C ./builddir install

singularity --version
