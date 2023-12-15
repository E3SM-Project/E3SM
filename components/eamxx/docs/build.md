# Building and Testing SCREAM Unit Tests

Follow these simple instructions to build and test SCREAM's standalone
configuration for yourself. Note that similar documentation is available on confluence (for E3SM team members) 
at https://acme-climate.atlassian.net/wiki/spaces/NGDNA/pages/1264386127/Running+SCREAM+Tests. 
This document makes use of the following paths:

+ `${RUN_ROOT_DIR}`: the root directory where SCREAM is built and run
+ `${SCREAM_SRC_DIR}`: the directory into which you've cloned the `scream` repo

SCREAM's configuration and build system is based on [CMake](https://cmake.org/).
CMake has been around a while and has gained a lot of traction in recent years,
especially in the HPC community. It has good [reference documentation](https://cmake.org/cmake/help/latest/index.html),
but it can be tricky to use if you've never encountered it. Ask a SCREAM team
member for help if you're stuck on something CMake-related.

If you see a `CMakeLists.txt` files or a file with a `.cmake` suffix, that's
just part of the build system. You might also see files with `CTest` as part of
their name. These files are related to [CTest](https://cmake.org/cmake/help/latest/manual/ctest.1.html),
CMake's testing tool.

## 1. Start From a Trustworthy Commit

First, make sure you've cloned the [SCREAM repo (including all submodules)](https://github.com/E3SM-Project/scream)
to `SCREAM_SRC_DIR` using the following command:

```
git clone --recurse-submodules https://github.com/E3SM-Project/scream
```

If you have already cloned the project and forgot to type `--recurse-submodules`,
you can change to `$SCREAM_SRC_DIR` and using the following command to initialize,
fetch and checkout all submodules:

```
git submodule update --init --recursive
```

If you're running a branch that's not `master`, check out this branch with

```
git checkout <branch>
```

## 2. Configure Your SCREAM Build

Change to your `$RUN_ROOT_DIR` directory and use CMake to configure your build.

If you're building SCREAM on one of our supported platforms, you can tell CMake
to use the appropriate machine file using the `-C` flag. Machine files are
located in `$SCREAM_SRC_DIR/components/eamxx/cmake/machine-files`. Take a look
and see whether your favorite machine has one.

For example, to configure SCREAM on the Quartz machine at LLNL:

```
cd $RUN_ROOT_DIR
cmake \
    -DCMAKE_CXX_COMPILER=$(which mpicxx) \
    -DCMAKE_BUILD_TYPE=Debug \
    -C ${SCREAM_SRC_DIR}/components/eamxx/cmake/machine-files/quartz.cmake \
    ${SCREAM_SRC_DIR}/components/eamxx
```

If you're building on a machine that doesn't have a ready-made machine file,
you can try configuring your build by manually passing options to CMake. This
usually looks something like the following:
```
cd $RUN_ROOT_DIR
cmake \
    -D CMAKE_BUILD_TYPE=Debug \
    -D Kokkos_ENABLE_DEBUG=TRUE \
    -D Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION=OFF \
    -D Kokkos_ENABLE_SERIAL=ON \
    -D Kokkos_ENABLE_OPENMP=ON \
    -D Kokkos_ENABLE_PROFILING=OFF \
    -D Kokkos_ENABLE_DEPRECATED_CODE=OFF \
    -D KOKKOS_ENABLE_ETI:BOOL=OFF \
    -D CMAKE_C_COMPILER=mpicc \
    -D CMAKE_CXX_COMPILER=mpicxx \
    -D CMAKE_Fortran_COMPILER=mpif90 \
    ${SCREAM_SRC_DIR}/components/eamxx
```

In either case, SCREAM requires MPI-savvy compilers, which can be specified
using the `CMAKE_xyz_COMPІLER` options.

Above, we've configured `Debug` builds to make it easier to find and fix errors.
For performance testing, you should configure a `Release` build and make use of
other options, depending on your architecture.

## 3. Build SCREAM

Now you can build SCREAM from that same directory:

```
make -j
```

The `-j` flag tells Make to use threads to compile in parallel. If you like, you
can set the number of threads by passing it as an argument to `-j` (e.g.
`make -j8`).

## 4. Run SCREAM's Tests

Before running the tests, generate a baseline file:

```
cd $RUN_ROOT_DIR
make baseline
```

The tests will run, automatically using the baseline file, which is located in
the CMake-configurable path `${SCREAM_TEST_DATA_DIR}`. By default, this path is
set to `data/` within your build directory (which is `$RUN_ROOT_DIR`, in
our case).

To run all of SCREAM's tests, make sure you're in `$RUN_ROOT_DIR` and type

```
ctest -VV
```

This runs everything and reports results in an extra-verbose (`-VV`) manner.

You can also run subsets of the SCREAM tests. For example, to run only the
P3 regression tests (again, from the `$RUN_ROOT_DIR` directory), use

```
ctest -R p3_regression
```

### Grouping Tests with Labels

We can create groupings of tests by using **labels**. For example, we have a
`driver` label that runs tests for SCREAM's standalone driver. You can see a
list of available labels by typing

```
ctest --print-labels
```

To see which tests are associated with a given label (e.g. `driver`), use

```
ctest -L driver -N
```

# SCREAM Test Suites

## The `p3_regression` Suite

`p3_regression` uses a baseline file to compare any new or altered
implementations with our P3 Fortran reference implementation. If you're working
on the C++/Kokkos implementation, you can invoke any new tests to the function
`Baseline::run_and_cmp` in
`${SCREAM_SRC_DIR}/components/eamxx/p3/tests/p3_run_and_cmp.cpp`.

If the reference Fortran implementation changes enough that a new baseline file
is required, make sure to let other SCREAM team members know, in order to
minimize disruptions.

