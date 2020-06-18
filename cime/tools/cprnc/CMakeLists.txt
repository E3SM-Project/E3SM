# Generate this with: $cimeroot/tools/configure --mpilib=mpi-serial --macros-format=CMake
#   You'll also need to source the .env_mach_specific.sh file before trying to build cprnc

include("${BLDROOT}/Macros.cmake")
set(CMAKE_C_COMPILER "${SCC}")
set(CMAKE_Fortran_COMPILER "${SFC}")

PROJECT(CPRNC C Fortran)
ENABLE_LANGUAGE(Fortran)
set(CMAKE_Fortran_FLAGS "${FFLAGS}")

message("HERE fortran flags are ${CMAKE_Fortran_FLAGS} FFLAGS are ${FFLAGS}")

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Find netcdf
SET(NetCDF_PATH ${NETCDF_PATH})
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{CIMEROOT}/src/externals/pio2/cmake)
find_package (NetCDF COMPONENTS Fortran REQUIRED)

ADD_CUSTOM_COMMAND(
  OUTPUT ${PROJECT_BINARY_DIR}/compare_vars_mod.F90
  COMMAND perl ${PROJECT_SOURCE_DIR}/../../src/externals/genf90/genf90.pl
    ${PROJECT_SOURCE_DIR}/compare_vars_mod.F90.in > ${PROJECT_BINARY_DIR}/compare_vars_mod.F90
  DEPENDS ${PROJECT_SOURCE_DIR}/compare_vars_mod.F90.in ${PROJECT_SOURCE_DIR}/../../src/externals/genf90/genf90.pl
)

execute_process(COMMAND ${NETCDF_PATH}/bin/nf-config --flibs
  RESULT_VARIABLE NF_CONFIG_RESULT
  OUTPUT_VARIABLE NF_CONFIG_OUTPUT)

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

if (NF_CONFIG_RESULT STREQUAL "0")
  separate_arguments(NF_CONFIG_LIST UNIX_COMMAND "${NF_CONFIG_OUTPUT}")
  target_link_libraries(cprnc ${NF_CONFIG_LIST})
else()
  target_link_libraries(cprnc ${NetCDF_Fortran_LIBRARIES} )
endif()
