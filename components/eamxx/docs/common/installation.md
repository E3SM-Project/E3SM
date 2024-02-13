# Installation

Follow these simple instructions to build and test EAMxx's standalone
configuration for yourself. This document makes use of the following paths:

+ `${RUN_ROOT_DIR}`: the root directory where EAMxx is built and run
+ `${EAMXX_SRC_DIR}`: the directory into which you've cloned the `scream` repo

EAMxx's configuration and build system is based on [CMake](https://cmake.org/).
CMake has been around a while and has gained a lot of traction in recent years,
especially in the HPC community. It has good [reference documentation](https://cmake.org/cmake/help/latest/index.html),
but it can be tricky to use if you've never encountered it. Ask a EAMxx team
member for help if you're stuck on something CMake-related.

If you see a `CMakeLists.txt` files or a file with a `.cmake` suffix, that's
just part of the build system. You might also see files with `CTest` as part of
their name. These files are related to [CTest](https://cmake.org/cmake/help/latest/manual/ctest.1.html),
CMake's testing tool.

## Prerequisites

First, make sure you're on one of the machines supported by EAMxx, or that you
have the following software installed:

* A working MPI installation (typically [MPICH]() or [Open-MPI]())
* [CMake](https://cmake.org) and [GNU Make](https://www.gnu.org/software/make/)
* A working set of C, C++, and Fortran compilers
* A recent version of [Git](https://git-scm.com/)
* A working installation of [NetCDF](https://www.unidata.ucar.edu/software/netcdf/),
  including both [C](https://github.com/Unidata/netcdf-c) and
  [Fortran](https://github.com/Unidata/netcdf-fortran) libraries.

## Setting Up Your Environment

## Configuring and Building Scream

### 1. Start From a Trustworthy Commit

First, make sure you've cloned the [EAMxx repo (including all submodules)](https://github.com/E3SM-Project/scream)
to `EAMXX_SRC_DIR` using the following command:

```
git clone --recurse-submodules https://github.com/E3SM-Project/scream
```

If you have already cloned the project and forgot to type `--recurse-submodules`,
you can change to `$EAMXX_SRC_DIR` and using the following command to initialize,
fetch and checkout all submodules:

```
git submodule update --init --recursive
```

If you're running a branch that's not `master`, check out this branch with

```
git checkout <branch>
```

### 2. Configure Your EAMxx Build

Change to your `$RUN_ROOT_DIR` directory and use CMake to configure your build.

If you're building SCREAM on one of our supported platforms, you can tell CMake
to use the appropriate machine file using the `-C` flag. Machine files are
located in `$EAMXX_SRC_DIR/components/eamxx/cmake/machine-files`. Take a look
and see whether your favorite machine has one.

For example, to configure SCREAM on the Quartz machine at LLNL:

```
cd $RUN_ROOT_DIR
cmake \
    -DCMAKE_CXX_COMPILER=$(which mpicxx) \
    -DCMAKE_BUILD_TYPE=Debug \
    -C ${EAMXX_SRC_DIR}/components/eamxx/cmake/machine-files/quartz.cmake \
    ${EAMXX_SRC_DIR}/components/eamxx
```

If you're building on a machine that doesn't have a ready-made machine file,
you can try configuring your build by manually passing options to CMake. This
usually looks something like the following, which configures EAMxx to compile
CPU code using Kokkos's OpenMP backend:
```
cd $RUN_ROOT_DIR
cmake \
    -D CMAKE_BUILD_TYPE=Debug \
    -D CMAKE_C_COMPILER=mpicc \
    -D CMAKE_CXX_COMPILER=mpicxx \
    -D CMAKE_Fortran_COMPILER=mpif90 \
    -D MPIEXEC_EXECUTABLE=`which mpiexec` \
    -D EKAT_MPI_NP_FLAG:STRING=-n \
    -D SCREAM_DYNAMICS_DYCORE=HOMME \
    -D SCREAM_DOUBLE_PRECISION:BOOL=ON \
    -D SCREAM_INPUT_ROOT:PATH=/path/to/scream-input \
    -D Kokkos_ENABLE_DEBUG=TRUE \
    -D Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION=OFF \
    -D Kokkos_ENABLE_SERIAL=ON \
    -D Kokkos_ENABLE_OPENMP=ON \
    -D Kokkos_ENABLE_LIBDL=OFF \
    -D Kokkos_ENABLE_PROFILING=OFF \
    -D Kokkos_ENABLE_DEPRECATED_CODE=OFF \
    -D KOKKOS_ENABLE_ETI:BOOL=OFF \
    -D NetCDF_C_PATHS=/path/to/netcdf-c-dir \
    -D NetCDF_Fortran_PATHS=/path/to/netcdf-f90-dir \
    -D PnetCDF_C_PATHS=/path/to/pnetcdf-dir \
    -D PnetCDF_Fortran_PATHS=/path/to/pnetcdf-f90-dir \
    ${EAMXX_SRC_DIR}/components/eamxx
```

In either case, EAMxx requires MPI-aware compilers. Let's examine these
options (only some of which are required on any given machine) to make sure we
know what they do:

* `CMAKE_BUILD_TYPE`: specifies whether you are building EAMxx in a
  developer-friendly configuration (`Debug`), for a production run (`Release`)
  or for performance profiling or some other specialized purpose. Typically,
  you'll set this option to `Debug` or `Release`.
* `CMAKE_{C,CXX,Fortran}_COMPILER`: the name of the command used to invoke an
  MPI-enabled C, C++, or Fortran compiler to build EAMxx
* `MPIEXEC_EXECUTABLE`: the name of the command used to run EAMxx using MPI,
  typically `mpiexec` or `mpirun`, but possibly different depending on your
  desired machine
* `EKAT_MPI_NP_FLAG`: the flag passed to `MPIEXEC_EXECUTABLE` that you use to
  specify the number of desired MPI processes. This is typically `-n` for
  `mpiexec` and `-np` for `mpirun`.
* `SCREAM_DYNAMICS_DYCORE`: specifies the dycore used for configuring EAMxx,
  which is `NONE` if you are not configuring EAMxx to run its dycore-related
  tests, or `HOMME` if you want to use HOMMExx
* `SCREAM_DOUBLE_PRECISION`: indicates whether EAMxx's `Real` type is a
  double-precision (`ON`) or single-precision (`OFF`) floating point type
* `SCREAM_INPUT_ROOT`: specifies the location of the top-level folder that
  stores input data files for EAMxx. This folder is populated with input files
  which are downloaded automatically during EAMxx's build process.
* The Kokkos-related build options (most of which begin with `Kokkos_`) are
  described [in the Kokkos Wiki](https://kokkos.github.io/kokkos-core-wiki/keywords.html)
* `NetCDF_C_PATHS`: specifies one or more folders in which the NetCDF C library
  and headers are installed. In the simplest configuration, the headers should
  be located in `${NetCDF_C_PATHS}/include` and the library should live in
  `${NetCDF_C_PATHS}/lib`.
* `NetCDF_Fortran_PATHS`: specifies one or more folders in which the NetCDF
  Fortran library and modules are installed. Analogous to `${NetCDF_C_PATHS}`,
  `.mod` files should be in `${NetCDF_Fortran_PATHS}/include`, and the library
  should be installed in `${NetCDF_Fortran_PATHS}/lib`.
* `PnetCDF_C_PATHS`: specifies one or more folders in which the pNetCDF C
  library and headers are installed, analogous to `NetCDF_C_PATHS`.
* `PnetCDF_Fortran_PATHS`: specifies one or more folders in which the pNetCDF
  Fortran library and modules are installed, analogous to
  `NetCDF_Fortran_PATHS`.

Above, we've configured `Debug` builds to make it easier to find and fix errors.
For performance testing, you should configure a `Release` build and make use of
other options, depending on your architecture.

### 3. Build SCREAM

Now you can build SCREAM from that same directory:

```
make -j
```

The `-j` flag tells Make to use threads to compile in parallel. If you like, you
can set the number of threads by passing it as an argument to `-j` (e.g.
`make -j8`).

## Running Tests

You can run EAMxx's tests to make sure your build works by following the
instructions [here](../developer/standalone_testing.md).
