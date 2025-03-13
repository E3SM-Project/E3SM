# Installation

You can build and run RDycore on the following platforms:

* Linux and Mac laptops and workstations
* [Frontier](https://www.olcf.ornl.gov/frontier/) (Oak Ridge National Laboratory)
* [Perlmutter](https://docs.nersc.gov/systems/perlmutter/) (NERSC)

## Required Software

To build RDycore, you need:

* [CMake v3.14+](https://cmake.org/)
* GNU Make
* reliable C, C++, and Fortran compilers
* a working MPI installation (like [OpenMPI](https://www.open-mpi.org/)
  or [MPICH](https://www.mpich.org/))
* [PETSc](https://petsc.org/release/), built with the following third-party
  libraries:
    * cgns
    * exodusii
    * fblaslapack
    * hdf5
    * libceed
    * metis
    * muparser
    * netcdf
    * parmetis
    * pnetcdf
    * zlib

You can obtain all of these freely on the Linux and Mac platforms. On Linux,
just use your favorite package manager. On a Mac, you can get the Clang C/C++
compiler by installing XCode, and then use a package manager like
[Homebrew](https://brew.sh/) or [MacPorts](https://www.macports.org/) to get
the rest.

### Which version of PETSc?

Check [our automated testing workflow](https://github.com/RDycore/RDycore/blob/main/.github/workflows/auto_test.yml#L24)
for the proper Git hash to use to build RDycore. The linked line specifies a
Docker image containing the "blessed" version of PETSc, which can be read as
follows:

```
coherellc/rdycore-petsc:fc288817-int32
```

* `coherellc` is the name of the DockerHub organization hosting the image
* `rdycore-petsc` is the the name of the Docker image
* `fc288817` is the Git hash within the [PETSc repository](https://gitlab.com/petsc/petsc)
  used to build RDycore
* `int32` (or `int64`) indicates whether the PETSc installation within the image
  uses 32-bit or 64-bit integers for the `PetscInt` data type.

See our [PETSc Dockerfile](https://github.com/RDycore/RDycore/blob/main/tools/Dockerfile.petsc#L50)
for an example of the `configure` command we use to build PETSc in our continous
integration environment.

## Clone the Repository

First, go get the [source code](https://github.com/RDycore/RDycore)
at GitHub:

=== "SSH"
    ```bash
    git clone git@github.com:RDycore/RDycore.git
    ```
=== "HTTPS"
    ```bash
    git clone https://github.com/RDycore/RDycore.git
    ```

This places an `RDycore` folder into your current path.

## Configure RDycore

RDycore uses CMake, and can be easily configured as long as
[PETSc is installed](https://petsc.org/release/install/) and [the `PETSC_DIR`
and `PETSC_ARCH` environment variables are set
properly](https://petsc.org/release/install/multibuild/#environmental-variables-petsc-dir-and-petsc-arch).
Usually all you need to do is change to your `RDycore` source directory and type

```bash
cmake -S . -B build
```

where `build` is the name of your build directory relative to the source
directory. If you want to install RDycore somewhere afterward, e.g. to be able
to configure E3SM to use it, you can set the prefix for the installation path
using the `CMAKE_INSTALL_PREFIX` parameter:

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/path/to/install
```

### Supported configuration options

CMake allows you to specify build options with the `-D` flag, as indicated in
Step 3 above. Here are the options supported by RDycore:

* **`CMAKE_INSTALL_PREFIX=/path/to/install`**: a path to which the RDycore library
  and driver are installed with `make install` (as in Step 7)
* **`CMAKE_BUILD_TYPE=Debug|Release`**: controls whether a build has debugging
  information or whether it is optimized
* **`CMAKE_VERBOSE_MAKEFILE=ON|OFF`**: if `ON`, displays compiler and linker
  output while building. Otherwise displays only the file being built.
* **`ENABLE_COVERAGE=ON|OFF`**: if `ON`, enables code coverage instrumentation.

Since RDycore gets most of its configuration information from PETSc, we don't
need to use most other CMake options.

### Considerations for Apple hardware

If you're on a Mac, make sure you have installed the XCode Command Line Tools.
If you have, these tools should be located in
`/Library/Developer/CommandLineTools/usr/bin/`, so add this directory to your
`PATH`.

## Build, Test, and Install RDycore

After you've configured RDycore, you can build it:

1. Change to your build directory (e.g. `cd build`)
2. Type `make -j` to build the library.
3. To run tests for the library (and the included drivers), type `make test`.
4. To install the model to the location (indicated by your `CMAKE_INSTALL_PREFIX`,
   if you specified it), type `make install`. By default, products are installed
   in the `include`, `lib`, `bin`, and `share` subdirectories of this prefix.

### Running Tests

RDycore uses [CTest](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html),
CMake's testing program, to run its tests. CTest is very fancy and allows us to
run tests selectively and in various ways, but all you need to do to run all the
tests for RDycore is to change to your build directory and type

```
make test
```

This runs every test defined in your build configuration and dumps the results
to `Testing/Temporary/LastTest.log`.

### Measuring Code Coverage

RDycore can use [gcov](https://gcc.gnu.org/onlinedocs/gcc/Gcov.html) or
[lcov](https://lcov.readthedocs.io/en/latest/index.html) to analyze code
coverage (the fraction of source code that is exercised by programs and tests)
with the GCC or Clang compilers.

To instrument the `rdycore` library and unit tests for code coverage analysis,
pass the `-DENABLE_COVERAGE=ON` flag to CMake when configuring your build. Then,
after building and running tests, type

```
make coverage
```

to generate a single report (`coverage.info`) containing all coverage
information. See the documentation for `gcov` and `lcov` (linked above) for
details on how to interpret thіs information.

### Checking for memory errors and leaks with Valgrind

If you're using a Linux system and have [Valgrind](https://valgrind.org/)
installed, you can run our tests using Valgrind's `memcheck` tool with

```
make memcheck
```

## Making code changes and rebuilding

Notice that you must build RDycore in a  **build tree**, separate from its source
trees. This is standard practice in CMake-based build systems, and it allows you
to build several different configurations without leaving generated and compiled
files all over your source directory. However, you might have to change the way
you work in order to be productive in this kind of environment.

When you make a code change, make sure you build from the build directory that
you created in step 1 above:

```bash
cd /path/to/RDycore/build
make -j
```

You can also run tests from this build directory with `make test`.

This is very different from how some people like to work. One method of making
this easier is to use an editor in a dedicated window, and have another window
open with a terminal, sitting in your `build` directory. If you're using a fancy
modern editor, it might have a CMake-based workflow that handles all of this for
you.

The build directory has a structure that mirrors the source directory, and you
can type `make` in any one of its subdirectories to do partial builds. In
practice, though, it's safest to always build from the top of the build tree.


## Preinstalled PETSc for RDycore on certain DOE machines

The RDycore team supports installation of the model at following DOE machines:

1. [Perlmutter](https://docs.nersc.gov/systems/perlmutter/) at NERSC
2. [Frontier](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html) at OLCF

First, run the following shell script to set PETSc-related environmental variables and load appropriate modules.

```bash
`source config/set_petsc_settings.sh --mach <machine_name> --config <configuration>`,
```

Multiple configurations of PETSc have been pre-installed on these supported machines
under RDycore's project directories. Information about the available PETSc configurations
can be obtained via `./config/set_petsc_settings.sh`.

The Perlmutter system has two types of compute nodes: CPU-only and CPU-GPU nodes, and
RDycore needs to be build separately for each type of compute node. The CPU-only nodes
have 128 cores (2 x 64-core AMD EPYC CPUs), while the CPU-GPU nodes have
1 x 64-core AMD EPYC CPU and 4 x NVIDIA A100. RDycore uses PETSc's and libCEED's support
of CUDA to run on Perlmutter GPUs.

Frontier has a single type of compute node that has 64-core AMD and 4x AMD MI250X GPUs.
Each GPU has 2 Graphics Compute Dies (GCDs) for a total of 8 GCDs per node. Of the 64-cores,
only 56 are allocatable cores instead of 64 cores. RDycore uses PETSc's and libCEED's
support of HIP to run on AMD GPUs.

### Example: Building and running RDycore on Perlmutter CPU nodes

```bash
cd /path/to/RDycore

# Set PETSc environment variables for Perlmutter CPU nodes
source config/set_petsc_settings.sh --mach pm-cpu --config 1

# Build RDycore
cmake -S . -B build-$PETSC_ARCH -DCMAKE_INSTALL_PREFIX=$PWD/build-$PETSC_ARCH
cd build-$PETSC_ARCH
make -j4 install

# Use an interactive job queue
salloc --nodes 1 --qos interactive --time 00:30:00 --constraint cpu \
--account=<project-id>

# Change to the directory containing tests
cd driver/tests/swe_roe

# Run on 4 MPI tasks on CPUs
srun -N 1 -n 4 ../../rdycore ex2b_ic_file.yaml -ceed /cpu/self -log_view
```

### Example: Building and running RDycore on Perlmutter GPU nodes

```bash
cd /path/to/RDycore

# Set PETSc environment variables for Perlmutter GPU nodes
source config/set_petsc_settings.sh --mach pm-gpu --config 1

# Build RDycore
cmake -S . -B build-$PETSC_ARCH -DCMAKE_INSTALL_PREFIX=$PWD/build-$PETSC_ARCH
cd build-$PETSC_ARCH
make -j4 install

# Use an interactive job queue
salloc --nodes 1 --qos interactive --time 00:30:00 --constraint gpu \
--gpus 4 --account=<project-id>_g

# Change to the directory containing tests
cd driver/tests/swe_roe

# Run on 4 GPUs using CUDA
srun -N 1 -n 4 -c 32 ../../rdycore ex2b_ic_file.yaml \
-ceed /gpu/cuda -dm_vec_type cuda -log_view -log_view_gpu_time
```

### Example: Building and running RDycore on Frontier

```bash
cd /path/to/RDycore

# Set PETSc environment variables for Frontier
source config/set_petsc_settings.sh --mach frontier --config 1

# Build RDycore
cmake -S . -B build-$PETSC_ARCH -DCMAKE_INSTALL_PREFIX=$PWD/build-$PETSC_ARCH
cd build-$PETSC_ARCH
make -j4 install

# Use an interactive job queue
salloc -N 1 -A <project-id> -t 0:30:00 -p batch

# Change to the directory containing tests
cd driver/tests/swe_roe

# Run on CPUs
srun -N 1 -n8 -c1 ../../rdycore ex2b_ic_file.yaml -ceed /cpu/self -log_view

# Run on 8 GPUs using HIP
srun -N 1 -n8 -c1 ../../rdycore ex2b_ic_file.yaml \
-ceed /gpu/hip -dm_vec_type hip -log_view -log_view_gpu_time
```
