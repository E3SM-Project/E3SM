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

## Dependencies
Get Kokkos:
```
cd ${KOKKOS_SRC_LOC}
git clone http://github.com/kokkos/kokkos
```
Configure it as a debug build and build it:
```
cd ${RUN_ROOT_DIR}/kokkos/build
rm -rf CMakeCache.txt CMakeFiles
cmake \
    -D CMAKE_INSTALL_PREFIX=${PWD}/../install \
    -D CMAKE_BUILD_TYPE=Debug \
    -D KOKKOS_ENABLE_DEBUG=ON \
    -D KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION=OFF \
    -D KOKKOS_ENABLE_SERIAL=ON \
    -D KOKKOS_ENABLE_OPENMP=ON \
    -D KOKKOS_ENABLE_PROFILING=OFF \
    -D KOKKOS_ENABLE_DEPRECATED_CODE=OFF \
    -D KOKKOS_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF \
    ../kokkos
make -j8 install
```
For performance testing, other options should be used depending on architecture.

## SCREAM
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

## Testing
From a trusted SHA1, generate a baseline file:
```
cd ${RUN_ROOT_DIR}/test
make baseline
```
The generated baseline file is in `${RUN_ROOT_DIR}/test/p3/tests`.

From any SHA1, run tests:
```
cd ${RUN_ROOT_DIR}/test
ctest -VV
```
The test `ctest -R p3_regression` uses the baseline file to compare any new or
altered implementations with the Fortran reference implementation's output in
the baseline file.

The intended use of the baseline capability is as follows:

1. Run `make baseline` from a trusted commit, such as master HEAD.

2. Develop code. If you're creating a new P3 implementation, add it to the
function `Baseline::run_and_cmp` in
`${SCREAM_SRC_LOC}/scream/components/scream/p3/tests/p3_run_and_cmp.cpp`.

3. The test `ctest -R p3_regression` will compare output against the baseline
file.

4. If the reference Fortran impl is changed such that output changes, carefully
decide whether that is desired. If it is, then in your commit, tell everyone
that a new baseline is needed.

It is not important to know where the baseline file is, but if it is of
interest, it is in the CMake-configurable `${SCREAM_TEST_DATA_DIR}`, which
defaults to `data/` in the build directory.
