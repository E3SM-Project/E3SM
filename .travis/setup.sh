#!/bin/bash

# Remove the pre-installed Go 1.11 from the environment
go version
for p in $(xargs -n1 -d: <<< $PATH); do
    echo $p | grep '/go' > /dev/null
    if [ $? -ne 0 ]; then
        TMPPATH="$p:$TMPPATH"
    fi
done
export PATH="$TMPPATH"
unset GOPATH
unset GOROOT

# Install Go 1.13.5
export VERSION=1.13.5 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz

export GOPATH="${HOME}/go"
export PATH="/usr/local/go/bin:${PATH}:${GOPATH}/bin"
go version

# Install Singularity
export VERSION=3.5.3 && # adjust this as necessary \
    wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz && \
    tar -xzf singularity-${VERSION}.tar.gz && \
    cd singularity
./mconfig && \
    make -C ./builddir && \
    sudo make -C ./builddir install

singularity --version
