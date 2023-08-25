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

## Standalone Build

CMake and OMEGA prefer an out-of-source build. This enables a user to build
and maintain multiple executables from the same source directory.
The standard practice is for a user to create a separate directory where
the build should take place and the commands below should be launched from
that directory.

To utilize a compiler name that's defined in the E3SM machine configuration
file (`config_machines.xml`), employ the `OMEGA_CIME_COMPILER` CMake variable,
as illustrated below:

```sh
>> cmake \
  -DOMEGA_CIME_COMPILER=crayclang \
  ${E3SM_HOME}/components/omega
```

To employ a specific compiler, users can utilize the `OMEGA_CXX_COMPILER`
CMake variable, providing CMake with the compiler's path.

```sh
>> cmake \
  -DOMEGA_CXX_COMPILER=CC \
  ${E3SM_HOME}/components/omega
```

`OMEGA_CXX_COMPILER` overrides `OMEGA_CIME_COMPILER`.

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
src
test
```

To build the Omega, execute the `make` command in the build directory.

Typical output will look something like:

```sh
>> make
Scanning dependencies of target yakl
[ 10%] Building Fortran object external/YAKL/CMakeFiles/yakl.dir/src/YAKL_gator_mod.F90.o
[ 20%] Building CXX object external/YAKL/CMakeFiles/yakl.dir/src/YAKL.cpp.o
[ 30%] Linking CXX static library libyakl.a
[ 30%] Built target yakl
[ 40%] Building CXX object src/CMakeFiles/OmegaLib.dir/base/MachEnv.cpp.o
[ 50%] Building CXX object src/CMakeFiles/OmegaLib.dir/ocn/OcnDummy.cpp.o
[ 60%] Linking CXX static library libOmegaLib.a
[ 60%] Built target OmegaLib
[ 70%] Building CXX object src/CMakeFiles/omega.exe.dir/drivers/DrvDummy.cpp.o
[ 80%] Linking CXX executable omega.exe
[ 80%] Built target omega.exe
[ 90%] Building CXX object test/CMakeFiles/testDataTypes.exe.dir/base/DataTypesTest.cpp.o
[100%] Linking CXX executable testDataTypes.exe
[100%] Built target testDataTypes.exe
```

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

To run the Omega test suite, execute the `ctest` command.

```sh
>> ctest
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
