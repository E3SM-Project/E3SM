# This cmake utility allows you to build cprnc just by adding
#
#   Include(BuildCprnc)
#   BuildCprnc()
#
# to your CMakeLists.txt.

function(BuildCprnc)

  set(BLDROOT ${PROJECT_BINARY_DIR}/externals/cprnc)
  file(WRITE ${BLDROOT}/Macros.cmake
"
set(SCC ${CMAKE_C_COMPILER})
set(SFC ${CMAKE_Fortran_COMPILER})
set(FFLAGS \"${CMAKE_Fortran_FLAGS}\")
set(NETCDF_PATH ${NetCDF_Fortran_PATHS})
"
)
  set(SRC_ROOT ${PROJECT_SOURCE_DIR}/../..)
  set(ENV{CIMEROOT} ${SRC_ROOT}/cime)
  add_subdirectory($ENV{CIMEROOT}/tools/cprnc ${BLDROOT})
endfunction()
