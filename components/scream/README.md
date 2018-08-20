# SCREAM standalone: How to set up basic testing.

## Dependencies
Get Kokkos:
```
git clone http://github.com/kokkos/kokkos
```
Configure it as a debug build and build it:
```
cd kokkos
cmake \
    -D CMAKE_INSTALL_PREFIX=/path/to/kokkos/install \
    -D CMAKE_BUILD_TYPE=Debug \
    -D KOKKOS_ENABLE_DEBUG=ON \
    -D KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION=OFF \
    -D KOKKOS_ENABLE_SERIAL=ON \
    -D KOKKOS_ENABLE_OPENMP=ON \
    -D KOKKOS_ENABLE_PROFILING=OFF \
    -D KOKKOS_ENABLE_DEPRECATED_CODE=OFF \
    -D KOKKOS_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF \
    ./
make -j8 install
```
For performance testing, other options should be used depending on architecture.

## SCREAM
Configure and build it:
```
cmake \
    -D Kokkos_DIR=/path/to/kokkos/install \
    -D SCREAM_ENABLE_FPE=FALSE \
    -D SCREAM_DOUBLE_PRECISION=FALSE \
    -D CMAKE_BUILD_TYPE=RelWithDebInfo \
    ~/climate/scream/components/scream
make -j8
```

## Testing
From a trusted SHA1, generate a baseline file:
```
make baseline
```
From any SHA1, run tests:
```
ctest -VV
```
