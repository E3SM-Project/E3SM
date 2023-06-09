include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

set(NetCDF_Fortran_PATH /usr/gdata/climdat/libs/netcdf-fortran/install/lassen/fortran CACHE STRING "")
set(BLAS_LIBRARIES /usr/gdata/climdat/libs/blas/libblas.a CACHE STRING "")
set(LAPACK_LIBRARIES /usr/gdata/climdat/libs/lapack/liblapack.a CACHE STRING "")

set(SCREAM_INPUT_ROOT "/usr/gdata/climdat/ccsm3data/inputdata/" CACHE STRING "")
