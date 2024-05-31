#!/bin/env bash

cmake \
    -S ../../../ \
    -B ../build_src \
    -DCMAKE_BUILD_TYPE='Release' \
    -DEAMXX_ENABLE_PYBIND='ON' \
    -DNetCDF_Fortran_PATH=$CONDA_PREFIX \
    -DNetCDF_C_PATH=$CONDA_PREFIX \
    -DCMAKE_CXX_FLAGS='-fvisibility-inlines-hidden -fmessage-length=0' \
    -DCMAKE_C_FLAGS='' \
    -DCMAKE_Fortran_FLAGS='' \
    -DCMAKE_CXX_COMPILER=mpic++ \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_Fortran_COMPILER=mpif90 \
    -DBUILD_SHARED_LIBS=ON \
    -DSCREAM_INPUT_ROOT="../build_src" \
    -DSCREAM_ENABLE_BASELINE_TESTS=OFF \
    -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX

cmake --build ../build_src/src/python -j${CPU_COUNT:-128}

find ../build_src -type f -name "*.so" | xargs cp -t ../libpyeamxx

python -m pip install -vvv ../
