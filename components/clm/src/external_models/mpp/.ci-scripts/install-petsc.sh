#!/bin/sh

git clone https://bitbucket.org/petsc/petsc.git

PETSC_GIT_HASH=a12052c5384dbd298e0348d2e27951245b7274e9

cd petsc

git checkout ${PETSC_GIT_HASH}

export PETSC_DIR=$PWD
export PETSC_ARCH=linux-gnu

./configure PETSC_ARCH=linux-gnu --with-debug=$DEBUG --with-shared-libraries=1 --download-hdf5 --download-metis --download-parmetis --with-c2html=0 --download-mpich

make all
make test

