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

First, make sure you've cloned the [SCREAM repo](https://github.com/E3SM-Project/scream)
to `SCREAM_SRC_DIR` (including all submodules). If you're running a branch that's not `master`, check out
this branch with

```
git checkout <branch>
```

## 2. Configure Your SCREAM Build

Change to your `RUN_ROOT_DIR` directory and use CMake to configure your build.
This usually looks something like the following:
```
cd $RUN_ROOT_DIR
cmake \
    -D CMAKE_BUILD_TYPE=Debug \
    -D KOKKOS_ENABLE_DEBUG=TRUE \
    -D KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION=OFF \
    -D KOKKOS_ENABLE_SERIAL=ON \
    -D KOKKOS_ENABLE_OPENMP=ON \
    -D KOKKOS_ENABLE_PROFILING=OFF \
    -D KOKKOS_ENABLE_DEPRECATED_CODE=OFF \
    -D KOKKOS_ENABLE_EXPLICIT_INSTANTIATION:BOOL=OFF \
    ${SCREAM_SRC_DIR}/components/scream
```

If you're building on your local laptop or workstation, make sure you have MPI
compilers installed, and tell SCREAM about them by inserting these options before
the last line of the above command:

```
    -D CMAKE_C_COMPILER=mpicc \
    -D CMAKE_CXX_COMPILER=mpicxx \
    -D CMAKE_Fortran_COMPILER=mpif90 \
```

Here, we've configured a `Debug` build to make it easier to find and fix errors.
For performance testing, you should configure a `Release` build and make use of
other options, depending on your architecture.

## 3. Build SCREAM

Now you can build SCREAM from that same directory:

```
make -j
```

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

To run all of SCREAM's tests, make sure you're in `RUN_ROOT_DIR` and type

```
ctest -VV
```

This runs everything and reports results in an extra-verbose (`-VV`) manner.

You can also run subsets of the SCREAM tests. For example, to run only the
P3 regression tests (again, from the `RUN_ROOT_DIR` directory), use

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
`${SCREAM_SRC_DIR}/components/scream/p3/tests/p3_run_and_cmp.cpp`.

If the reference Fortran implementation changes enough that a new baseline file
is required, make sure to let other SCREAM team members know, in order to
minimize disruptions.

