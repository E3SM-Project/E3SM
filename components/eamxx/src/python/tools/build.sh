#!/bin/env bash

cd components/eamxx/src/python

cmake \
    -S ../../ \
    -B build_src \
    -DCMAKE_BUILD_TYPE='Release' \
    -DEAMXX_ENABLE_PYBIND='ON' \
    -DNETCDF_FORTRAN_PATH=$PREFIX \
    -DNETCDF_C_PATH=$PREFIX \
    -DCMAKE_CXX_FLAGS='-fvisibility-inlines-hidden -fmessage-length=0 -Wno-use-after-free -Wno-unused-variable -Wno-maybe-uninitialized' \
    -DCMAKE_C_FLAGS='' \
    -DCMAKE_Fortran_FLAGS='-Wno-maybe-uninitialized -Wno-unused-dummy-argument' \
    -DCMAKE_CXX_COMPILER=mpic++ \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_Fortran_COMPILER=mpif90 \
    -DBUILD_SHARED_LIBS=ON \
    -DSCREAM_INPUT_ROOT="build_src" \
    -DSCREAM_ENABLE_BASELINE_TESTS=OFF \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DPYTHON_EXECUTABLE=$PYTHON

cmake --build build_src/src/python -j${CPU_COUNT:-128}

find build_src -type f -name "*.so" | xargs cp -t libpyeamxx/

$PYTHON -m build -w -n -x

pip install dist/*.whl -vv
