(omega-user-build)=

# CMake-based Omega Build

Omega's build system is constructed using the CMake build tool,
offering a strong foundation for managing the building process.

The Omega build system has two modes: standalone and E3SM component.

At the start of the Omega build, the system reads the E3SM machine file
(`config_machines.xml`) for both standalone and E3SM component builds.
Following this, it configures CMake variables and environment variables
based on the computing system where the build is taking place, as well as
user input from the CMake command-line.

For the Omega build system to function, a Python interpreter is necessary.
The minimum version of CMake is 3.21.

## Standalone Build

CMake and OMEGA prefer an out-of-source build. This enables a user to build
and maintain multiple executables from the same source directory.
The standard practice is for a user to create a separate directory where
the build should take place and the commands below should be launched from
that directory.

To utilize a compiler name that's defined in the E3SM machine configuration
file (`config_machines.xml`), employ the `OMEGA_CIME_COMPILER` CMake variable,
as illustrated below. In some cases, you may also want to add `OMEGA_CIME_MACHINE`
to specify which system you intend to use. The values of `OMEGA_CIME_COMPILER`
and `OMEGA_CIME_MACHINE` are defined in
"${E3SM}/cime\_config/machines/config\_machines.xml".
OMEGA requires some external libraries. Many of these are built automatically
from the E3SM distribution. However, the METIS, ParMETIS, and, optionally,
GKlib libraries must be built separately and the path must be supplied during
the cmake invocation as shown below. If a `OMEGA_METIS_ROOT` is not supplied,
it is assumed that both METIS and ParMETIS are installed in the same
`OMEGA_PARMETIS_ROOT` location. If a `OMEGA_GKLIB_ROOT` is not supplied, it is
assumed that both GKLIB and METIS are installed in the same `OMEGA_METIS_ROOT`
location. All three libraries should be static libraries.

```sh
>> cmake \
  -DOMEGA_CIME_COMPILER=nvidiagpu \
  -DOMEGA_CIME_MACHINE=pm-gpu \
  -DOMEGA_PARMETIS_ROOT=/path/to/parmetis \
  ${E3SM_HOME}/components/omega
```

Once the command completes successfully, several scripts will be created
in the build directory.

* omega\_env.sh   : load specific modules and set env. variables read from CIME
* omega\_build.sh : run `make` command after sourcing `omega_env.sh`
* omega\_run.sh   : run `./src/omega.exe` after sourcing `omega_env.sh`
* omega\_ctest.sh : run `ctest` after sourcing `omega_env.sh`

Run omega\_build.sh in the build directory to build Omega.

To employ a specific compiler, users can utilize the `OMEGA_CXX_COMPILER`
CMake variable, providing CMake with the compiler's path.

```sh
>> cmake \
  -DOMEGA_CXX_COMPILER=CC \
  ${E3SM_HOME}/components/omega
```

`OMEGA_CXX_COMPILER` overrides `OMEGA_CIME_COMPILER`.

To build Omega for GPU, add OMEGA\_ARCH to one of "CUDA" or "HIP".

```sh
>> cmake \
  -DOMEGA_CIME_COMPILER=nvidiagpu \
  -DOMEGA_CIME_MACHINE=pm-gpu \
  -DOMEGA_ARCH=CUDA \
  ${E3SM_HOME}/components/omega
```

To enable the ctest-based unittest option, include the `OMEGA_BUILD_TEST`
option as illustrated below.

```sh
>> cmake \
  -DOMEGA_BUILD_TEST=ON \
  ${E3SM_HOME}/components/omega
```

Once the `cmake` command succeeds, the directory where the command was
executed will contain the following files and directories.

```sh
>> ls -1 ${WORKDIR}
CMakeCache.txt
CMakeFiles
cmake_install.cmake
CTestTestfile.cmake
external
Makefile
omega_build.sh
omega_ctest.sh
omega_env.sh
omega_run.sh
src
test
```

To build the Omega, execute the `./omega_build.sh` command in the build directory.

If the build succeeds, the Omega library and executable are created in the
`src` sub-directory of the build directory.

To specify the output directory for saving the generated output files,
a user can include the `OMEGA_INSTALL_PREFIX` variable in the `cmake` command,
as shown below:

```sh
cmake \
  ... \
  -DOMEGA_INSTALL_PREFIX="/path/to/output/directory" \
  ${E3SM_HOME}/components/omega
  ...
```

Afterwards, a user can use the `make install` command to copy the Omega library
and executable to the `lib` and `bin` sub-directories under the directory
specified by `DOMEGA_INSTALL_PREFIX`.

The `./omega_ctest.sh` command runs Omega unit tests. To run the tests, MPI
parallel job launcher such as SLURM srun should be available. You may first
get allocation of an interactive computing node or use batch system.
In addition, you must either copy (or soft link) a mesh file into the same
directory as the unit test executables or specify in a configuration file
(not implemented yet) the path to that mesh input file. For detailed instructions
on how to obtain appropiate test meshes see {ref}`omega-dev-quick-start-getting-meshes`.

For example, if you are on an interactive computing node, you can run
Omega unit test by running `omega_ctest.sh` as shown below.

```sh
>> ./omega_ctest.sh
Test project <cmake working directory>
    Start 1: DATA_TYPES_TEST
1/1 Test #1: DATA_TYPES_TEST ..................   Passed    0.03 sec

100% tests passed, 0 tests failed out of 1

Total Test time (real) =   0.04 sec
```

## E3SM Component Build

In the E3SM component build, the compilation is triggered by
the E3SM build system.

### Manual Preparation for E3SM Component Build

Although all the build configurations in the E3SM component build
should ideally be managed through the E3SM build configuration,
the current Omega build is not yet integrated into the E3SM build.
Therefore, users need to make specific modifications to enable
the Omega build within the E3SM build process.

Step 1: Modify `${E3SM_ROOT}/components/CMakeLists.txt`

Add the lines indicated as "added" in the `CMakeLists.txt` file.

```bash
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_util.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_mpas_model.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_omega.cmake) # <= added
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_eamxx.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_model.cmake)
...
set(BUILDCONF ${CASEROOT}/Buildconf)

build_mpas_models()
build_omega() # <= added
```

Step 2: Create `${E3SM_ROOT}/components/cmake/build_omega.cmake`

Create the `cmake` file and copy the following content into it.

```cmake
function(build_omega)

  # Set CIME source path relative to components
  set(CIMESRC_PATH "../cime/src")

  add_subdirectory("omega")

endfunction(build_omega)
```

Step 3: Create an E3SM Case

At this stage, a user can create any E3SM case as usual without
requiring any specific configuration for Omega.

NOTE: In this version, the compiled library file for Omega,
"libOmegaLib.a", is not linked to the E3SM executable. You can
find the Omega library file at `${E3SM_BLD_DIRECTORY}/cmake-bld/omega/src`.
