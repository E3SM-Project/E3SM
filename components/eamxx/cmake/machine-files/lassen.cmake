include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

set(NetCDF_PATH /g/g20/foucar1/packages/netcdf CACHE STRING "")
set(LAPACK_LIBRARIES /usr/lib64/liblapack.so CACHE STRING "")
set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

set(SCREAM_INPUT_ROOT "/usr/gdata/climdat/ccsm3data/inputdata/" CACHE STRING "")
