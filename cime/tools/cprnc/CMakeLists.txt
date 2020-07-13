# Generate this with: $cimeroot/tools/configure --mpilib=mpi-serial --macros-format=CMake
#   You'll also need to source the .env_mach_specific.sh file before trying to build cprnc

include("${BLDROOT}/Macros.cmake")
set(CMAKE_C_COMPILER "${SCC}")
set(CMAKE_Fortran_COMPILER "${SFC}")

project(CPRNC C Fortran)
enable_language(Fortran)
set(CMAKE_Fortran_FLAGS "${FFLAGS}")

message("HERE fortran flags are ${CMAKE_Fortran_FLAGS} FFLAGS are ${FFLAGS}")

cmake_minimum_required(VERSION 2.8)

# Find netcdf
set(NetCDF_PATH ${NETCDF_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{CIMEROOT}/src/externals/pio2/cmake)
find_package (NetCDF COMPONENTS Fortran REQUIRED)

# generate compare_vars_mod.F90
add_custom_command(
  OUTPUT ${PROJECT_BINARY_DIR}/compare_vars_mod.F90
  COMMAND perl ${PROJECT_SOURCE_DIR}/../../src/externals/genf90/genf90.pl
    ${PROJECT_SOURCE_DIR}/compare_vars_mod.F90.in > ${PROJECT_BINARY_DIR}/compare_vars_mod.F90
  DEPENDS ${PROJECT_SOURCE_DIR}/compare_vars_mod.F90.in ${PROJECT_SOURCE_DIR}/../../src/externals/genf90/genf90.pl
)

# Set up includes
include_directories(
  ${NetCDF_Fortran_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
)

#
# Set up lib dependencies, relying on nf-config if possible or NetCDF_Fortran_LIBRARIES
# as a fallback. We want the executable to find libraries via RPATH so that cprnc is
# less-dependent on the current environment (since cprnc is built with a serial netcdf,
# it's likely that a parallel case will have different netcdf modules loaded when it
# comes time to execute cprnc).
#
execute_process(COMMAND ${NETCDF_PATH}/bin/nf-config --flibs
  RESULT_VARIABLE NF_CONFIG_RESULT
  OUTPUT_VARIABLE NF_CONFIG_OUTPUT)

if (NF_CONFIG_RESULT STREQUAL "0")
  separate_arguments(NF_LIB_LIST UNIX_COMMAND "${NF_CONFIG_OUTPUT}")
else()
  set(NF_LIB_LIST ${NetCDF_Fortran_LIBRARIES})
endif()

message("lib list is: ${NF_LIB_LIST}")

foreach(NF_LIB IN LISTS NF_LIB_LIST)
  if (NF_LIB MATCHES "-l")
    continue()
  elseif (NF_LIB MATCHES "-L/")
    string(REGEX REPLACE "^-L" "" NF_LIB_DIR "${NF_LIB}")
    list(APPEND NF_LIB_DIRS ${NF_LIB_DIR})
  else()
    get_filename_component(NF_LIB_DIR ${NF_LIB} DIRECTORY)
    list(APPEND NF_LIB_DIRS ${NF_LIB_DIR})
  endif()
endforeach()

message("lib dirs are: ${NF_LIB_DIRS}")

set(CMAKE_BUILD_RPATH ${NF_LIB_DIRS})

# Add targets
set (CPRNC_SRCS
  ${PROJECT_BINARY_DIR}/compare_vars_mod.F90
  filestruct.F90
  utils.F90
  prec.F90
  cprnc.F90
)

add_executable(cprnc ${CPRNC_SRCS})
add_dependencies(cprnc COMPARE_VARS)

target_link_libraries(cprnc ${NF_LIB_LIST})
