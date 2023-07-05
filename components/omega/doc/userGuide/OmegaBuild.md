(omega-user-build)=

# CMake-based Omega Build

Omega's build system is built upon the CMake build tool, which provides
a robust foundation for managing the build process.

The Omega build system supports two modes: standalone and E3SM component.

## Standalone Build

CMake and OMEGA prefer an out-of-source build. This enables a user to build
and maintain multiple executables from the same source directory.
The standard practice is for a user to create a separate directory where
the build should take place and the commands below should be launched from
that directory.

To perform a standalone build, you need to execute the cmake command with
the required CMake and Omega parameters. For example, you can specify the
CC parameter as the C++ compiler and enable the ctest option to guide the
Omega build process. The following example demonstrates building
a standalone Omega with a test suite.

```sh
>> cmake \
  -DOMEGA_BUILD_TEST=ON \
  -DOMEGA_CXX_COMPILER=CC \
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

To build the Omega, execute the `make` command in the directory.

```sh
>> make
Scanning dependencies of target yakl
[ 11%] Building Fortran object external/YAKL/CMakeFiles/yakl.dir/src/YAKL_gator_mod.F90.o
[ 22%] Building CXX object external/YAKL/CMakeFiles/yakl.dir/src/YAKL.cpp.o
[ 33%] Linking CXX static library libyakl.a
[ 33%] Built target yakl
[ 44%] Building CXX object src/CMakeFiles/OmegaLib.dir/ocn/ocndummy.cpp.o
[ 55%] Linking CXX static library libOmegaLib.a
[ 55%] Built target OmegaLib
[ 66%] Building CXX object src/CMakeFiles/omega.exe.dir/drivers/drvdummy.cpp.o
[ 77%] Linking CXX executable omega.exe
[ 77%] Built target omega.exe
[ 88%] Building CXX object test/CMakeFiles/OMEGA_TEST.dir/testdummy.cpp.o
[100%] Linking CXX executable OMEGA_TEST
[100%] Built target OMEGA_TEST
```

If the build succeeds, the Omega library and executable are created in the
`src` sub-directory of the build directory.

To specify the output directory for saving the generated output files,
you can include the `OMEGA_INSTALL_PREFIX` variable in the `cmake` command,
as shown below:

```sh
cmake \
  ... \
  -DOMEGA_INSTALL_PREFIX="/path/to/output/directory" \
  ${E3SM_HOME}/components/omega
  ...
```

Afterwards, you can use the `make install` command to copy the Omega library
and executable to the `lib` and `bin` sub-directories under the directory
specified by `DOMEGA_INSTALL_PREFIX`.

To run the Omega test suite, execute the `ctest` command.

```sh
>> ctest
Test project <cmake working directory>
    Start 1: OMEGA_TEST
1/1 Test #1: OMEGA_TEST .......................   Passed    0.01 sec

100% tests passed, 0 tests failed out of 1

Total Test time (real) =   0.01 sec
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

At this stage, you can create any E3SM case as usual without
requiring any specific configuration for Omega.

NOTE: In this version, the compiled library file for Omega,
"libOmegaLib.a", is not linked to the E3SM executable. You can
find the Omega library file at `${E3SM_BLD_DIRECTORY}/cmake-bld/omega/src`.
