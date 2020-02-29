# Generate this with: $cimeroot/tools/configure --mpilib=mpi-serial --macros-format=CMake
#   You'll also need to source the .env_mach_specific.sh file before trying to build cprnc

include("$ENV{BLDROOT}/Macros.cmake")
set(CMAKE_Fortran_Compiler "${SFC}")
set(CMAKE_Fortran_FLAGS "${FFLAGS}")

PROJECT(CPRNC C Fortran)

ENABLE_LANGUAGE(Fortran)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Find netcdf
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{CIMEROOT}/src/externals/pio2/cmake)
find_package (NetCDF COMPONENTS Fortran REQUIRED)

ADD_CUSTOM_COMMAND(
  OUTPUT ${PROJECT_BINARY_DIR}/compare_vars_mod.F90
  COMMAND perl ${PROJECT_SOURCE_DIR}/../../src/externals/genf90/genf90.pl
    ${PROJECT_SOURCE_DIR}/compare_vars_mod.F90.in > ${PROJECT_BINARY_DIR}/compare_vars_mod.F90
  DEPENDS ${PROJECT_SOURCE_DIR}/compare_vars_mod.F90.in ${PROJECT_SOURCE_DIR}/../../src/externals/genf90/genf90.pl
)

INCLUDE_DIRECTORIES(
  ${NetCDF_Fortran_INCLUDE_DIRS}
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

TARGET_LINK_LIBRARIES(cprnc ${NetCDF_Fortran_LIBRARIES})
