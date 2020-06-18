
cmake_minimum_required (VERSION 2.8.12)

# Adjust CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

option(CMAKE_INSTALL_PREFIX "Sets installation prefix [/usr/local].")
option(BUILD_SHARED_LIBS "Builds shared libraries [ON]." ON)
#option(WITH_NETCDF "Enables support for netcdf [OFF]." ON)
#option(TPL_NETCDF_LIBRARIES "List of absolute paths to netcdf link libraries [].")
#option(TPL_NETCDF_INCLUDE_DIRS "List of absolute paths to netcdf include directories [].")

# For now, we disable shared libs on Macs.
if (APPLE)
  set(BUILD_SHARED_LIBS OFF)
endif()

if (NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX /usr/local)
endif()

#set(PETSC_ARCH $ENV{PETSC_ARCH})
#set(PETSC_DIR $ENV{PETSC_DIR})
#include($ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/initial_cache_file.cmake)

# Make sure compilers are set. This must be done before enabling languages.
#if (NOT CMAKE_C_COMPILER)
#  if (NOT $ENV{CC} STREQUAL "")
#    set(CMAKE_C_COMPILER $ENV{CC})
#  else()
#    set(CMAKE_C_COMPILER cc)
#  endif()
#endif()
#if (NOT CMAKE_C_FLAGS)
#  set(CMAKE_C_FLAGS $ENV{CFLAGS})
#endif()
#if (NOT CMAKE_CXX_COMPILER)
#  if (NOT $ENV{CXX} STREQUAL "")
#    set(CMAKE_CXX_COMPILER $ENV{CXX})
#  else()
#    set(CMAKE_CXX_COMPILER c++)
#  endif()
#endif()
#if (NOT CMAKE_CXX_FLAGS)
#  set(CMAKE_CXX_FLAGS $ENV{CXX_FLAGS})
#endif()
if (NOT CMAKE_Fortran_COMPILER)
  if (NOT $ENV{FC} STREQUAL "")
    set(CMAKE_Fortran_COMPILER $ENV{FC})
  else()
    set(CMAKE_Fortran_COMPILER gfortran)
  endif()
endif()
if (NOT CMAKE_Fortran_FLAGS)
  set(CMAKE_Fortran_FLAGS $ENV{FCFLAGS})
endif()
#enable_language(C)
#enable_language(CXX)
enable_language(Fortran)

# We declare the project here.
project (EMI)

#message("-- C compiler is ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})")
#message("-- CXX compiler is ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID})")
message("-- Fortran compiler is ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})")

if (BUILD_SHARED_LIBS)
  message("-- EMI will be built as a shared library.")
else()
  message("-- EMI will be built as a static library.")
endif()

set (EMI_VERSION_MAJOR 0)
set (EMI_VERSION_MINOR 1)
set (EMI_VERSION_PATCH 0)
set (EMI_VERSION "${EMI_MAJOR_VERSION}.${EMI_MINOR_VERSION}.${EMI_PATCH_VERSION}")

# General C compiler flags.
#if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -pedantic-errors -Wextra -Werror-implicit-function-declaration")
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")
#
#  if (BUILD_SHARED_LIBS)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fPIC")
#  endif()
#
#  if (LINUX EQUAL 1)
#    # Counter some of GCC's more recent stinginess on Linux.
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_POSIX_C_SOURCE=200809L")# -D_BSD_SOURCE")
#  endif()
#
#elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang")
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -pedantic-errors -Wextra -Werror-implicit-function-declaration -fno-builtin")
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -Wno-unused-function")
#elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")
#endif()
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")

# Fortran compiler flags.
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -W -Wall -std=gnu -pedantic -ffree-line-length-0 -Wno-unused-variable -Wno-unused-parameter -DCPRGNU -fprofile-arcs -ftest-coverage")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCPRINTEL")
endif()


# Figure out the system type.
set(BETR_HAVE_BOOL 1) # All reasonable C99 compilers have this now.
if (APPLE EQUAL 1)
  set(SYS_FLAGS "-DAPPLE=1")
  set(DYLIB_SUFFIX "dylib")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
else ()
  if (LINUX EQUAL 1)
    set(SYS_FLAGS "-DLINUX=1")
    set(DYLIB_SUFFIX "so")
  else()
    if (WIN32 EQUAL 1)
      set(BETR_HAVE_BOOL 0) # MS doesn't have reasonable C compilers.
      set(SYS_FLAGS "-DWINDOWS=1")
      set(DYLIB_SUFFIX "dll")
    endif()
  endif ()
endif ()

# Here we make sure CMake-installed binaries use the correct runpath, and
# that the path is not stripped during installation.
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


#find_package(PETSc)
#if (NOT PETSC_FOUND)
#  message(FATAL_ERROR "PETSc was not found.")
#endif()
#include_directories(${PETSC_INCLUDES})
# Include PETSc in the rpath.
#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:${PETSC_DIR}/${PETSC_ARCH}/lib")

#find_package(NetCDF COMPONENTS F90 REQUIRED)
#if (NOT NETCDF_FOUND)
#  message(FATAL_ERROR "NetCDF was not found.")
#endif()
#include_directories(${NETCDF_INCLUDE_DIR})
#include_directories(${NETCDF_F90_INCLUDE_DIR})


# Other third-party libraries.
#add_subdirectory(3rd-party)

# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/include")
link_directories("${PROJECT_BINARY_DIR}/lib")

# Unit testing.
enable_testing()

# Source code itself.
include_directories("${PROJECT_SOURCE_DIR}")
add_subdirectory(elm_stub)
add_subdirectory(constants)
add_subdirectory(emi_data_types)
add_subdirectory(emi_data_definition)
add_subdirectory(emi_data_exchange)
add_subdirectory(em)
add_subdirectory(emi)
add_subdirectory(demo)

add_test(
NAME regression_tests
COMMAND regression_tests.py --backtrace --executable local/bin/standalone_EMI
WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/regression_tests
)



