PROJECT(CPRNC C Fortran)

# Generate this with: $cimeroot/tools/configure --mpilib=mpi-serial --macros-format=CMake
#   You'll also need to source the .env_mach_specific.sh file before trying to build cprnc
include("${PROJECT_BINARY_DIR}/Macros.cmake")
set(CMAKE_Fortran_Compiler "${SFC}")
set(CMAKE_Fortran_FLAGS "${FFLAGS}")

if (NOT NETCDF_PATH)
  message(FATAL_ERROR "Requires NETCDF_PATH either from Macros.cmake or argument")
endif()

ENABLE_LANGUAGE(Fortran)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{CIMEROOT}/src/externals/pio2/cmake)

ADD_CUSTOM_COMMAND(
  OUTPUT ${PROJECT_BINARY_DIR}/compare_vars_mod.F90
  COMMAND perl ${PROJECT_SOURCE_DIR}/../../src/externals/genf90/genf90.pl
    ${PROJECT_SOURCE_DIR}/compare_vars_mod.F90.in > ${PROJECT_BINARY_DIR}/compare_vars_mod.F90
)

ADD_CUSTOM_TARGET(COMPARE_VARS DEPENDS ${PROJECT_BINARY_DIR}/compare_vars_mod.F90)

INCLUDE_DIRECTORIES(
  ${NETCDF_PATH}/include
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
)

SET (CPRNC_SRCS
  ${PROJECT_BINARY_DIR}/compare_vars_mod.F90
  filestruct.F90
  utils.F90
  prec.F90
  cprnc.F90
)

ADD_EXECUTABLE(cprnc ${CPRNC_SRCS})
ADD_DEPENDENCIES(cprnc COMPARE_VARS)

# SLIBS comes from Macros.cmake
separate_arguments(SLIBS_LIST UNIX_COMMAND "${SLIBS}")
TARGET_LINK_LIBRARIES(cprnc ${SLIBS_LIST})
