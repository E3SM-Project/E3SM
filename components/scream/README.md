# SCREAM standalone: How to set up basic testing.

## Create directory structure
```
#${RUN_ROOT_DIR} is the root directory where code is built and run
#${KOKKOS_SRC_LOC} is the directory where kokkos is cloned into
#${SCREAM_SRC_LOC} is the directory where scream is cloned into
mkdir -p ${RUN_ROOT_DIR}/kokkos/build   #where kokkos is built
mkdir -p ${RUN_ROOT_DIR}/kokkos/install #where kokkos is installed
mkdir -p ${RUN_ROOT_DIR}/test/          #where scream is installed
                                        #and tests are run
```

## Create dependencies
Get Kokkos:
```
cd ${KOKKOS_SRC_LOC}
git clone http://github.com/kokkos/kokkos
#The next step avoids bug where kokkos_generated_settings.cmake
#can get put in unexpected directory locations:
git checkout develop
```
Configure it as a debug build and build it:
```
cd ${RUN_ROOT_DIR}/kokkos/build
rm -rf CMakeCache.txt CMakeFiles
cmake \
    -D CMAKE_INSTALL_PREFIX=${RUN_ROOT_DIR}/kokkos/install \
    -D CMAKE_BUILD_TYPE=Debug \
    -D KOKKOS_ENABLE_DEBUG=ON \
    -D KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION=OFF \
    -D KOKKOS_ENABLE_SERIAL=ON \
    -D KOKKOS_ENABLE_OPENMP=ON \
    -D KOKKOS_ENABLE_PROFILING=OFF \
    -D KOKKOS_ENABLE_DEPRECATED_CODE=OFF \
    -D KOKKOS_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF \
    ${KOKKOS_SRC_LOC}
make -j8 install
```
For performance testing, other options should be used depending on architecture.

## Build SCREAM
Configure and build it:
```
cd ${RUN_ROOT_DIR}/test
cmake \
    -D Kokkos_DIR=${RUN_ROOT_DIR}/kokkos/install \
    -D SCREAM_ENABLE_FPE=FALSE \
    -D SCREAM_DOUBLE_PRECISION=FALSE \
    -D CMAKE_BUILD_TYPE=RelWithDebInfo \
    ${SCREAM_SRC_LOC}/scream/components/scream
make -j8
```

## Run tests
After building SCREAM from a trusted SHA1 (such as master HEAD), generate
a baseline file:
```
cd ${RUN_ROOT_DIR}/test
make baseline
```
It is not necessary to know where the baseline file is, but if it is of
interest, it is in the CMake-configurable `${SCREAM_TEST_DATA_DIR}`, which
defaults to `data/` in the build directory.

Make your desired code modifications, then rebuild SCREAM from this new
code. At this point you can run all tests with:
```
cd ${RUN_ROOT_DIR}/test
ctest -VV
```
To just run the p3_regression test, for example, execute:
```
cd ${RUN_ROOT_DIR}/test/
ctest -R p3_regression
```

## p3_regression test details

p3_regression uses the baseline file to compare any new or altered
implementations with the Fortran reference implementation's output in
the baseline file. If you've created a new implementation (e.g. a
Kokkos version of P3 or a version testing out a new performance paradigm),
add the new implementation to the function `Baseline::run_and_cmp` in
`${SCREAM_SRC_LOC}/scream/components/scream/p3/tests/p3_run_and_cmp.cpp`.

If the reference Fortran impl is changed such that output changes, carefully
decide whether that is desired. If it is, then in your commit, tell everyone
that a new baseline is needed.

